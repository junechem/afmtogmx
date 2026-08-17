"""Read a ``pycryoff-ff/1`` JSON document into the same object a ``.off`` produces.

pycryoff is a Python reimplementation of CRYOFF's force-matching engine. Its own ``.off`` output
is a *report* — human-readable, free to change layout — so it publishes the fitted force field as
a versioned JSON document instead. This module maps that document onto exactly the attributes
:class:`~afmtogmx.core.gen_md.ReadOFF` builds from a CRYOFF ``.off``, so every backend
(``off.gmx``, ``off.openmm``, ``off.write_report``) works on it unchanged.

Three things arrive here that a ``.off`` cannot carry, and they are the reason this path exists:

* **Per-atom charges.** A force-matching fit determines charge *products* (``q_i q_j``), never
  individual charges, so ``ReadOFF`` initialises every charge to 0.0 and waits for a ``.charges``
  file. Forgetting that step silently writes a neutral-everywhere topology. The JSON carries the
  resolved per-atom charges when the fit produced them, and ``None`` when it did not — so the
  absence is explicit and this reader can warn instead of writing zeros.
* **Parameters are name-keyed.** No positional convention to agree on: the reader looks up ``V``
  and ``delta`` by name, rather than relying on ``NCO`` storing ``(V, m, delta)`` while ``COS``
  stores ``(delta, V, m)``.
* **pycryoff's extensions** — explicit polarizability (``[POL]``), combination rules (``[COMB]``),
  out-of-plane virtual sites — which CRYOFF has no syntax for at all. Polarizability is surfaced
  on ``off.polarization`` rather than being dropped, because a model fitted *with* induction is
  not the same model without it.

Units are unchanged from the ``.off`` path: kcal/mol, Angstrom, degrees. Conversion to kJ/nm
happens in the writers (``topology.py``, ``xml_generation.py``), as before.
"""
import json
import warnings

from . import functions

#: The schema major this reader implements. A document from a newer major is refused rather than
#: read best-effort — a silently mis-read force field is far worse than a failed load.
SUPPORTED_SCHEMA = "pycryoff-ff/1"

#: JSON ``card`` -> the key ``gen_empty_bonded()`` uses. Only the coupled dihedral differs.
_CARD_KEYS = {"BON": "BON", "ANG": "ANG", "DIH": "DIH", "BD3": "BD3", "CDI": "CDI",
              "CDIH": "CDI"}


def _check_schema(doc):
    schema = doc.get("schema", "")
    if schema != SUPPORTED_SCHEMA:
        raise ValueError(
            f"unsupported force-field schema {schema!r}; this build of afmtogmx reads "
            f"{SUPPORTED_SCHEMA!r}. Regenerate the document with a matching pycryoff, or "
            f"upgrade afmtogmx.")
    units = doc.get("units", {})
    # The writers assume the .off's native units. A document in kJ/nm would produce a force field
    # wrong by a factor of ~4184 with no other symptom, so check rather than trust.
    expected = {"length": "angstrom", "energy": "kcal/mol", "angle": "degree"}
    wrong = {k: units.get(k) for k, v in expected.items() if units.get(k) != v}
    if wrong:
        raise ValueError(
            f"force-field document is not in afmtogmx's native units: expected {expected}, "
            f"got {wrong}. Refusing rather than converting silently.")


def _param_tuple(term):
    """A term's parameters as the float tuple the bonded dicts key on.

    Order comes from the document, which lists them in the functional form's canonical (CRYOFF)
    order — the same order the ``.off`` parser produces, so both paths key identically and a
    consumer cannot tell them apart.
    """
    return tuple(float(p["value"]) for p in term["params"])


def _build_bonded(doc):
    """``{molname: {ATO, BON, ANG, BD3, DIH, CDI, EXC}}``, matching the ``.off`` parser's shape."""
    bonded = {}
    for mol in doc.get("molecules", []):
        entry = functions.gen_empty_bonded()
        for atom in mol.get("atoms", []):
            index = int(atom["index"])
            vdw, cou = atom["vdw_type"], atom["cou_type"]
            if atom.get("pseudo"):
                # Net-force / torque pseudo-atoms, normalized to upper case. Every filter in this
                # package tests `!= 'NETF'` exactly, and CRYOFF's .off echo happens to upper-case
                # them so that has always worked. A deck may spell them 'NetF'/'Torq', which the
                # JSON preserves — leaving that casing here would quietly carry two fitting
                # targets through into the charge table and the residue definitions as if they
                # were atoms.
                vdw = cou = atom["pseudo"]
            entry["ATO"]["All"][index] = (vdw, cou)
            vsite = atom.get("virtual_site")
            if vsite:
                # The .off parser stores the raw token list it split off the atom line. Keep a
                # structured equivalent instead: the rule is (weight, parent) pairs plus a kind,
                # and re-flattening it to tokens would only invite a second parser.
                entry["ATO"]["Virtual"][(index, atom["vdw_type"], atom["cou_type"])] = {
                    "kind": vsite.get("kind", "average"),
                    "rule": [(float(w), int(i)) for w, i in vsite.get("rule", [])],
                }

        for term in mol.get("bonded", []):
            card = _CARD_KEYS.get(term.get("card", "").upper())
            form = term.get("form", "").upper()
            if card is None or form not in entry[card]:
                warnings.warn(f"{mol['name']}: skipping unsupported bonded term "
                              f"{term.get('card')}/{term.get('form')}")
                continue
            # Param-major, like the .off path: the parameter tuple is the key and the value is
            # every atom group sharing it. Terms that fitted to identical values therefore merge,
            # which is what the GROMACS writers expect (one #define per distinct parameter set).
            entry[card][form].setdefault(_param_tuple(term), []).extend(
                [int(i) for i in group] for group in term.get("groups", []))

        entry["EXC"] = [[int(i) for i in pair] for pair in mol.get("exclusions", [])]
        bonded[mol["name"]] = entry
    return bonded


def _build_nonbonded(doc):
    """``{(typeA, typeB): {InteractionType: [params, ...]}}``, pair key sorted, ``COU*`` folded."""
    nonbonded = {}
    for term in doc.get("nonbonded", []):
        pair = term.get("pair")
        if not pair:
            continue
        key = tuple(sorted(pair))
        itype = term.get("raw_name", term["form"]).upper()
        if "COU" in itype:
            itype = "COU"                      # matches gen_md._gen_nonbonded
        nonbonded.setdefault(key, {}).setdefault(itype, []).append(
            [float(p["value"]) for p in term["params"]])
    return nonbonded


def _build_charges(doc, bonded):
    """``{molname: {atom_name: charge}}``, from the document rather than defaulted to zero."""
    per_atom = (doc.get("charges") or {}).get("per_atom")
    if per_atom is None:
        warnings.warn(
            "this force field carries no per-atom charges — the fit determined charge *products* "
            "and nothing resolved them into individual charges. All charges are 0.0; set them "
            "with off.load_charges_from_file() or off.charges before writing any topology.")
        per_atom = {}
    charges = {}
    for mol, entry in bonded.items():
        charges[mol] = {
            cou: float(per_atom.get(cou, 0.0))
            for _, (vdw, cou) in entry["ATO"]["All"].items()
            if vdw not in ("NETF", "TORQ")
        }
    return charges


def load_document(path):
    """Read and validate a schema document from ``path``. Returns the raw dict."""
    with open(path, "r") as handle:
        doc = json.load(handle)
    _check_schema(doc)
    return doc


def populate(obj, path):
    """Fill ``obj`` — a bare :class:`ReadOFF` — from the JSON document at ``path``.

    Sets exactly the attributes the ``.off`` parser sets, then leaves ``_finalize`` to build the
    charges, residues and backends the same way both paths do.
    """
    doc = load_document(path)
    obj.off_loc = path
    obj.document = doc
    obj._ff_bonded = {}
    obj.sections = {"ff_input": "", "intra_potential": "", "inter_potential": "",
                    "molecular_definition": "", "table_potential": ""}
    obj.bonded = _build_bonded(doc)
    obj.nonbonded = _build_nonbonded(doc)
    #: pycryoff extensions with no ``.off`` equivalent. Surfaced rather than dropped: a model
    #: fitted with explicit induction is not the same model without it, and a consumer that
    #: writes a plain point-charge topology from it is producing a different force field.
    obj.polarization = doc.get("polarization")
    obj.combinations = doc.get("combinations") or []
    obj.provenance = doc.get("provenance") or {}
    obj.fit = doc.get("fit")
    if obj.polarization is not None:
        warnings.warn(
            "this force field was fitted with explicit mutual polarization ([POL]); its "
            "electrostatics is permanent charges PLUS induced dipoles. Writing it as a fixed "
            "point-charge topology reproduces only the permanent part. See off.polarization.")
    return _build_charges(doc, obj.bonded)
