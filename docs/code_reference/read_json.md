# `read_json.py` — reading a pycryoff force field

Builds a `ReadOFF` from a `pycryoff-ff/1` JSON document instead of a CRYOFF `.off`.

[pycryoff](https://github.com/junechem/pycryoff) is a Python reimplementation of CRYOFF's
force-matching engine. Its own `.off` output is a human-readable *report* whose layout is not a
contract, so it publishes the fitted force field as a versioned JSON document alongside it. That
document is what this module reads.

Entry point is `ReadOFF.from_json(json_loc)` (`gen_md.py`); everything here is what it calls.

```python
import afmtogmx
off = afmtogmx.ReadOFF.from_json("inter.json")
off.openmm.gen_xml(output="forcefield.xml")
```

The resulting object is **indistinguishable** from one parsed out of a `.off` — same `bonded`,
`nonbonded`, `charges`, `residues`, same native units (kcal/mol, Angstrom, degrees) — so both
backends and `write_report` work on it unchanged. Unit conversion to kJ/nm still happens only in
the writers.

## Why this exists rather than a pycryoff-emitted `.off`

The `.off` path takes connectivity from the echoed deck and numbers from the `Intra-`/
`Inter-Potential` sections, and joins the two **by position** through the module-global counter
`functions.total_bonded_added` — no name check, no validation, and a mismatch assigns the wrong
parameters silently. It also has no keyword for anything CRYOFF lacked, so pycryoff's extensions
(`[POL]`, `[COMB]`, `[CHG]`, out-of-plane virtual sites) would be dropped without a word. A
polarizable model flattened to point charges is a *different* force field.

## Three things the JSON carries that a `.off` cannot

1. **Per-atom charges.** A force-matching fit determines charge *products* (`q_i q_j`), never
   individual charges, so `ReadOFF.__init__` sets every charge to `0.0` and waits for a
   `.charges` file. Forget that step and the topology is silently neutral everywhere. Here
   `charges.per_atom` supplies them — or is explicitly `null`, and the reader **warns**.
2. **Name-keyed parameters.** No positional convention to agree on: `NCO` stores `(V, m, delta)`
   while `COS` stores `(delta, V, m)`, and the document names each one.
3. **pycryoff's extensions**, surfaced on the object rather than dropped: `off.polarization`,
   `off.combinations`, plus `off.provenance` and `off.fit`.

## Module contents

### Constants
- `SUPPORTED_SCHEMA = "pycryoff-ff/1"` — the schema major this build reads.
- `_CARD_KEYS` — JSON `card` → the key `gen_empty_bonded()` uses (`CDIH` folds to `CDI`).

### Functions

`_check_schema(doc)`
- Raises `ValueError` if `schema` is not `SUPPORTED_SCHEMA`, or if `units` are not
  `{length: angstrom, energy: kcal/mol, angle: degree}`. Both are refusals, not conversions: a
  document in kJ/nm is wrong by ~4184 with no other symptom, and best-effort parsing of an unknown
  schema yields a silently wrong force field.

`_param_tuple(term)`
- A term's parameters as the float tuple the bonded dicts key on, in the document's order — the
  functional form's canonical (CRYOFF) order, identical to what the `.off` parser produces.

`_build_bonded(doc)` → `{molname: {ATO, BON, ANG, BD3, DIH, CDI, EXC}}`
- Starts from `functions.gen_empty_bonded()` so every key a consumer expects is present.
- `ATO['All'][atnum] = (vdw_type, cou_type)`, **int keys throughout** (the `.off` path mixes `int`
  for virtual atoms with `str` for normal ones; this path does not inherit that).
- **Normalizes `NETF`/`TORQ` to upper case.** Every filter in this package compares against
  `'NETF'` exactly, which works on a `.off` only because CRYOFF's echo happens to upper-case them.
  A deck may spell them `NetF`/`Torq` and the JSON preserves that — leaving it would carry two
  fitting targets into the charge table and the residue definitions as if they were atoms.
- `ATO['Virtual'][(atnum, at1, at2)]` is a dict `{'kind': 'average'|'oop', 'rule': [(w, parent)]}`
  — structured, rather than the `.off` path's raw token list.
- Bonded terms are inverted term-major → **param-major**: the parameter tuple keys a list of the
  atom groups sharing it, so terms fitted to identical values merge (one `#define` per distinct
  parameter set, which is what the GROMACS writers expect). Unsupported card/form combinations
  warn and are skipped.

`_build_nonbonded(doc)` → `{(typeA, typeB): {ITYPE: [params, ...]}}`
- Pair key `tuple(sorted(...))`; `ITYPE` is `raw_name` (so `EXPINTER`, `STRC` survive) except that
  anything containing `COU` folds to `'COU'`, matching `gen_md._gen_nonbonded`.

`_build_charges(doc, bonded)` → `{molname: {atom_name: charge}}`
- From `charges.per_atom`, NETF/TORQ excluded. Warns and falls back to `0.0` when the document
  carries none.

`load_document(path)` → validated raw `dict`.

`populate(obj, path)` → the charges dict
- Fills a bare `ReadOFF` with `off_loc`, `document`, `sections` (empty — there are no `.off` text
  sections), `bonded`, `nonbonded`, `polarization`, `combinations`, `provenance`, `fit`; returns
  the charges for `_finalize` to install. Warns when `polarization` is present, because writing a
  polarizable model as a fixed point-charge topology reproduces only its permanent part.

## Tests

`test/test_read_json.py`, against two fixtures vendored in `test/sample_json/`:

- `butanol_intra.json` — a real pycryoff atomic-force fit of a 1-butanol QM/MM cluster: 8 bond
  types, 13 angle types, 4 NCO dihedrals, 77 nonbonded pairs, two molecules.
- `extensions_pol.json` — a synthetic deck reaching `[POL]`, `[COMB]`, an `OOP` virtual site and
  declared per-atom charges.
