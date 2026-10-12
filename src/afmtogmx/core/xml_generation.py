"""XML section builders for OpenMM force-field files.

Module-level helpers that turn ``ReadOFF`` data structures (``bonded``,
``nonbonded``, ``charges``) into OpenMM ``<ForceField>`` XML sections.
Used by :class:`afmtogmx.core.openmm_backend.OpenMMBackend`.

Atom types are *qualified* by molecule name ("UNK_C1") to keep them
globally unique even when two ``.off`` molecules happen to reuse the
same raw type label.
"""
import math
import warnings
from collections import defaultdict

try:
    from openmm.app import Element as _OpenMMElement
    _OPENMM_AVAILABLE = True
except ImportError:
    _OPENMM_AVAILABLE = False


# --------------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------------

def _atom_map(ato_all):
    """Normalise ``ATO['All']`` to ``{int_id: (atname, attype)}``."""
    return {int(k): v for k, v in ato_all.items()}


def _qualify(mol, attype):
    """Return ``"<mol>_<attype>"`` — globally unique atom type name."""
    return f'{mol}_{attype}'


def _is_virtual(bonded, mol, attype):
    """True if ``attype`` belongs to a virtual site in ``mol``."""
    for (atnum, atname, vtype) in bonded[mol].get('ATO', {}).get('Virtual', {}):
        if vtype == attype:
            return True
    return False


def _element_symbol(attype):
    """Best-effort element symbol from an atom-type name.

    Atom-type names encode the element as a prefix (``OW`` -> oxygen,
    ``HW`` -> hydrogen). Two-letter elements such as sodium (``NA``) and
    chlorine (``CL``) only resolve correctly if both leading characters are
    considered, so this tries the first two characters as an element symbol
    first and falls back to the first character. Returns the canonical OpenMM
    symbol (e.g. ``'Na'``, ``'Cl'``, ``'O'``) on a match, otherwise the first
    character of the type name.
    """
    if _OPENMM_AVAILABLE:
        for n in (2, 1):
            if len(attype) >= n:
                try:
                    return _OpenMMElement.getBySymbol(attype[:n]).symbol
                except Exception:
                    continue
    return attype[0]


def _get_mass(element_symbol):
    if not _OPENMM_AVAILABLE:
        return 0.0
    try:
        el = _OpenMMElement.getBySymbol(element_symbol)
        return float(el.mass.value_in_unit(el.mass.unit))
    except Exception:
        return 0.0


def _unique_atom_names(atom_map):
    """Per-residue unique names like ``C0, C1, H0, H1`` (element + counter)."""
    counters, result = {}, {}
    for atid in sorted(atom_map.keys()):
        atname, attype = atom_map[atid]
        if attype in ('NETF', 'TORQ'):
            continue
        element = _element_symbol(attype)
        n = counters.get(element, 0)
        counters[element] = n + 1
        result[atid] = f'{element}{n}'
    return result


def _type_sort_key(name):
    """Order type names numerically when they are numbers, alphabetically otherwise."""
    return (0, int(name), '') if name.isdigit() else (1, 0, name)


def build_type_names(atom_types, numeric=False):
    """``{qualified_type: the name its ``<Type>`` carries}``.

    Normally the two are the same and this is the identity, which is what keeps the XML
    readable: a type is called ``UNK_C0``.

    ``numeric=True`` renames them ``"1"``, ``"2"``, ... and it is not a style choice.
    OpenMM's ``AmoebaMultipoleGenerator`` reads its parameters with
    ``int(atom.attrib['type'])`` and assigns them with
    ``atomTypes = [int(data.atomType[atom]) for atom in data.atoms]`` -- so the moment an
    ``<AmoebaMultipoleForce>`` appears anywhere in the file, EVERY atom type name in the
    force field has to parse as an integer or ``createSystem`` dies with a ``ValueError``
    on a name it never explains. Tinker's own AMOEBA ffxml files are numbered for exactly
    this reason.

    The descriptive name is not lost: it stays on the type's ``class``, and every other
    section refers to atoms by ``class`` rather than by type, so only ``<AtomTypes>`` and
    the ``<Residue>`` atom lines ever see the number.
    """
    if not numeric:
        return {qualified: qualified for _, _, qualified in atom_types}
    return {qualified: str(i + 1) for i, (_, _, qualified) in enumerate(atom_types)}


def _matrix(n, val=0.0):
    return [[val] * n for _ in range(n)]


def _matrix_str(matrix):
    return '\t'.join(str(v) for row in matrix for v in row)


# --------------------------------------------------------------------------
# Data collection — exposed for the orchestrator
# --------------------------------------------------------------------------

def collect_atom_types(bonded, mol_names):
    """Return list of ``(mol, raw_type, qualified_type)`` for every type used.

    Order follows the appearance of atoms in each molecule's ATO section.
    NETF/TORQ atoms are skipped.
    """
    seen, result = set(), []
    for mol in mol_names:
        for atnum, (atname, attype) in bonded[mol]['ATO']['All'].items():
            if attype in ('NETF', 'TORQ'):
                continue
            qualified = _qualify(mol, attype)
            if qualified not in seen:
                seen.add(qualified)
                result.append((mol, attype, qualified))
    return result


def build_type_to_charge(bonded, charges, atom_types):
    """Return ``{qualified_type: charge}`` joining ATO types with the charges dict.

    A charge belongs to the **Coulomb** type, which is ``attype`` — the third column of a CRYOFF
    atom line ``<label> <VDWtype> [COUtype]`` — and that is how both constructors key
    ``self.charges``: the ``.off`` path seeds it from ``pair[1]`` and
    :func:`read_json._build_charges` from the same field.

    This used to look the charge up by ``atname``, the **vdW** type, via a ``name_lookup`` that
    mapped each Coulomb type back to the first vdW type seen with it. On a deck where the two
    columns are equal — every deck before phenol — that is a no-op, which is why it survived. On
    a split-typed deck it silently returns 0.0 for every Coulomb type that is not also the name
    of a vdW type, and the resulting XML is complete, well-formed, and has no electrostatics on
    most of the molecule. Phenol's 5 vdW / 9 Coulomb typing kept charges only on ``O0``, ``C1``
    and ``H0`` — the three names common to both namespaces — and zeroed the six ortho/meta/para
    types.

    There is no caller-side workaround: three distinct charges (``C2_o``, ``C2_m``, ``C2_p``)
    cannot be stored under the single vdW key ``C2``.
    """
    return {qualified: charges.get(mol, {}).get(raw, 0.0)
            for mol, raw, qualified in atom_types}


def collect_nonbonded(nonbonded, atom_types, bonded=None):
    """Walk ``nonbonded`` and split entries by interaction type.

    Returns
    -------
    exp_entries  : list of ``(qual1, qual2, A, alpha)`` (kcal·Å units)
    str_entries  : list of ``(qual1, qual2, P1, P2, P3)`` (kcal·Å units)
    srd_by_power : dict ``{power: [(qual1, qual2, P1, r0)]}`` (kcal·Å units)
    cpn_entries  : list of ``(qual1, qual2, [11 params])`` (kcal·Å units), charge
                   penetration; **orientation-sensitive**, see :func:`gen_cpn_force`

    Raw .off type names (``OW``) are expanded to all matching qualified
    types (``H2OQM_OW``, ``H2OMM_OW``); shared parameters get applied to
    every qualified pair, matching the .off file's intent.

    Two namespaces meet here and they are not the same one. A CRYOFF atom line is
    ``<label> <VDWtype> [COUtype]``: the nonbonded cards (EXP/SRD/POW/...) are keyed on
    the **vdW** type, while ``atom_types`` — and therefore every OpenMM ``<Type>`` — is
    built from the **Coulomb** type, because that is what carries a per-atom charge.
    Most force fields never declare the third column, the two names coincide, and the
    distinction is invisible.  When they differ (say five ring carbons sharing one
    repulsion type but split into ortho/meta/para charges) matching them by name finds
    nothing: ``raw_to_qualified.get('C2')`` misses, the pair is dropped, and since the
    lookup tables in the section builders are zero-filled the result is an XML with
    ``A = 0`` for those pairs -- no exception, no warning, and no repulsion.

    So pass ``bonded`` to map each qualified type back to the vdW type of the atoms that
    use it.  Without it this falls back to assuming the two namespaces coincide, which is
    correct for every force field that omits the Coulomb column and wrong silently for
    the ones that do not.
    """
    raw_to_qualified = defaultdict(list)
    if bonded is None:
        for mol, raw, qualified in atom_types:
            raw_to_qualified[raw].append(qualified)
    else:
        seen = set()
        for mol, cou, qualified in atom_types:
            for _atid, (vdw, atom_cou) in bonded[mol]['ATO']['All'].items():
                if atom_cou != cou or vdw in ('NETF', 'TORQ'):
                    continue
                if (vdw, qualified) not in seen:      # one vdW type serves many Coulomb types
                    seen.add((vdw, qualified))
                    raw_to_qualified[vdw].append(qualified)

    exp_entries, str_entries, srd_by_power, cpn_entries = [], [], {}, []
    dropped = set()

    def _add_srd(q1, q2, P1, power, r0):
        srd_by_power.setdefault(power, []).append((q1, q2, P1, r0))

    for (raw1, raw2), params in nonbonded.items():
        # A pair carrying only COU is a charge-product term; the custom forces build nothing
        # from it and <NonbondedForce> carries the charges instead. Butanol's MM embedding
        # shell is entirely COU, so an unmatched name there is normal and must not warn --
        # only an unmatched name on a pair that *would* have produced a force is a defect.
        if any(t not in ('COU',) for t in params):
            for raw in (raw1, raw2):
                if raw not in raw_to_qualified:
                    dropped.add(raw)
        for q1 in raw_to_qualified.get(raw1, []):
            for q2 in raw_to_qualified.get(raw2, []):
                for itype, param_sets in params.items():
                    for pset in param_sets:
                        if itype == 'EXP':
                            exp_entries.append((q1, q2, float(pset[0]), float(pset[1])))
                        elif itype in ('STR', 'STRC'):
                            str_entries.append((q1, q2, float(pset[0]), float(pset[1]), float(pset[2])))
                        elif itype == 'SRD':
                            _add_srd(q1, q2, float(pset[0]), float(pset[1]), float(pset[2]))
                        elif itype == 'POW':
                            _add_srd(q1, q2, float(pset[0]), float(pset[1]), 0.0)
                        elif itype == 'CPN':
                            cpn_entries.append((q1, q2, [float(v) for v in pset]))
                        elif itype == 'BUC':
                            p1, p2, p3 = float(pset[0]), float(pset[1]), float(pset[2])
                            exp_entries.append((q1, q2, p1, p3))
                            _add_srd(q1, q2, p2, -6.0, 0.0)

    if dropped:
        warnings.warn(
            f"{len(dropped)} vdW type(s) named by the nonbonded cards match no atom type and "
            f"their interactions were discarded: {', '.join(sorted(dropped)[:8])}"
            f"{', ...' if len(dropped) > 8 else ''}. The lookup tables are zero-filled, so the "
            f"affected pairs would get zero repulsion and zero dispersion. This is what a vdW / "
            f"Coulomb type-name mismatch looks like: pass `bonded` so the two can be told apart.",
            stacklevel=2)

    return exp_entries, str_entries, srd_by_power, cpn_entries


# --------------------------------------------------------------------------
# Section builders
# --------------------------------------------------------------------------

def gen_atomtypes(bonded, atom_types, type_names=None):
    """``<AtomTypes>`` with one ``<Type>`` per qualified atom type.

    ``type_names`` (see :func:`build_type_names`) supplies the ``name``; the ``class`` is
    always the qualified type, so a numbered file still says what each type is.
    """
    names = type_names or build_type_names(atom_types)
    lines = ['<AtomTypes>']
    for mol, raw, qualified in atom_types:
        name = names[qualified]
        if _is_virtual(bonded, mol, raw):
            lines.append(f'<Type name="{name}" class="{qualified}" mass="0.0"/>')
        else:
            element = _element_symbol(raw)
            mass = _get_mass(element)
            lines.append(f'<Type name="{name}" class="{qualified}" element="{element}" mass="{mass}"/>')
    lines.append('</AtomTypes>')
    return '\n'.join(lines)


def gen_residues(bonded, mol_names, molname_translations, type_names=None):
    """``<Residues>`` listing atoms, bonds, and virtual sites for every molecule.

    A residue's atoms must name a ``<Type>``, not a class -- OpenMM binds one type per atom
    -- so this is the second and last place ``type_names`` is consulted.
    """
    names = type_names or {}
    lines = ['<Residues>']
    for mol in mol_names:
        resname  = molname_translations.get(mol, mol)
        atom_map = _atom_map(bonded[mol]['ATO']['All'])
        unique   = _unique_atom_names(atom_map)

        lines.append(f'<Residue name="{resname}">')

        for atid in sorted(atom_map.keys()):
            if atid not in unique:
                continue
            atname, attype = atom_map[atid]
            qualified = _qualify(mol, attype)
            lines.append(f'<Atom name="{unique[atid]}" '
                         f'type="{names.get(qualified, qualified)}"/>')

        written_bonds = set()
        for bond_type, pairs_dict in bonded[mol].get('BON', {}).items():
            for params, atom_pairs in pairs_dict.items():
                for at1_id, at2_id in atom_pairs:
                    if at1_id not in unique or at2_id not in unique:
                        continue
                    bond_key = tuple(sorted([at1_id, at2_id]))
                    if bond_key in written_bonds:
                        continue
                    written_bonds.add(bond_key)
                    lines.append(f'<Bond atomName1="{unique[at1_id]}" atomName2="{unique[at2_id]}"/>')

        for (atnum, atname, attype), definition in \
                bonded[mol].get('ATO', {}).get('Virtual', {}).items():
            atid = int(atnum)
            if atid not in unique:
                continue
            vsite_xml = _virtual_site_xml(unique[atid], definition, unique)
            if vsite_xml:
                lines.append(vsite_xml)
            else:
                # Never drop this quietly. Without a <VirtualSite> element OpenMM builds the
                # site as a free particle: it carries its charge, it is never repositioned,
                # and nothing downstream complains -- the model is simply not the model.
                raise ValueError(
                    f"{mol}: atom {atid} ({atname}) is a virtual site but no OpenMM "
                    f"<VirtualSite> could be written for it from {definition!r}. "
                    "Supported: 2- and 3-parent linear averages.")

        lines.append('</Residue>')
    lines.append('</Residues>')
    return '\n'.join(lines)


def _virtual_site_rule(definition):
    """Normalise either Virtual-site representation to ``[(weight, parent_id), ...]``.

    ``ATO['Virtual']`` carries two different shapes depending on which reader built it:

    * the ``.off`` text parser (``functions.parse_section``) stores the raw token tuple it
      split off the atom line, e.g. ``('3:', '0.6', '1', '+', '0.2', '2', '+', '0.2', '3')``;
    * the JSON reader (``read_json``) stores the structured form
      ``{'kind': 'average', 'rule': [(0.6, 1), (0.2, 2), (0.2, 3)]}``.

    This function used to understand only the first, so a force field loaded
    ``from_json`` silently lost every ``<VirtualSite>`` line — the topology still built,
    the M site was simply an unconstrained free particle. Handle both.
    """
    if isinstance(definition, dict):
        kind = definition.get('kind', 'average')
        if kind != 'average':
            # An out-of-plane site is NOT a linear average of its parents. Emitting
            # <VirtualSite type="average3"> for one would place it in the parents' plane
            # and look perfectly well-formed while being the wrong geometry.
            raise ValueError(f"virtual-site kind {kind!r} is not a linear average")
        return [(float(w), int(i)) for w, i in definition.get('rule', [])]
    parts     = [p for p in definition if p != '+']  # Strip '+' separators
    num_atoms = int(parts[0].rstrip(':'))
    rule, i = [], 1
    while i < len(parts) and len(rule) < num_atoms:
        rule.append((float(parts[i]), int(parts[i + 1])))
        i += 2
    return rule


def _virtual_site_xml(site_name, definition, unique):
    """Parse a virtual-site definition → ``<VirtualSite>`` element."""
    try:
        rule = _virtual_site_rule(definition)
        num_atoms = len(rule)
        weights   = [w for w, _ in rule]
        atom_ids  = [i for _, i in rule]
        ref_names = []
        for aid in atom_ids:
            if aid not in unique:
                return None
            ref_names.append(unique[aid])
        if num_atoms == 3:
            return (f'<VirtualSite type="average3" siteName="{site_name}" '
                    f'atomName1="{ref_names[0]}" atomName2="{ref_names[1]}" '
                    f'atomName3="{ref_names[2]}" weight1="{weights[0]}" '
                    f'weight2="{weights[1]}" weight3="{weights[2]}"/>')
        if num_atoms == 2:
            return (f'<VirtualSite type="average2" siteName="{site_name}" '
                    f'atomName1="{ref_names[0]}" atomName2="{ref_names[1]}" '
                    f'weight1="{weights[0]}" weight2="{weights[1]}"/>')
    except (ValueError, IndexError, TypeError, KeyError):
        return None


def gen_nonbonded_force(atom_types, type_to_charge, charges_elsewhere=False):
    """``<NonbondedForce>`` carrying point charges (sigma=epsilon=0).

    ``charges_elsewhere=True`` zeroes every charge here because an
    ``<AmoebaMultipoleForce>`` is carrying the permanent electrostatics instead. The force
    is still written: it is inert, but it is where the periodic method and the cutoff are
    declared, and ``afm_openmm.prepare_afm_system`` reads both off it. Leaving the charges
    in as well would double every Coulomb interaction in the system.
    """
    lines = ['<NonbondedForce coulomb14scale="1.0" lj14scale="1.0">']
    for mol, raw, qualified in atom_types:
        charge = 0.0 if charges_elsewhere else type_to_charge.get(qualified, 0.0)
        lines.append(f'<Atom class="{qualified}" charge="{charge}" sigma="0.0" epsilon="0.0"/>')
    lines.append('</NonbondedForce>')
    return '\n'.join(lines)


def gen_bond_force(bonded, mol_names):
    """Unified ``<CustomBondForce>`` for HAR (k3=k4=0) and QUA bonds.

    Returns ``''`` if no bonds are present.
    """
    entries, written = [], set()

    for mol in mol_names:
        if 'BON' not in bonded[mol]:
            continue
        atom_map = _atom_map(bonded[mol]['ATO']['All'])

        for (r0, k2, k3, k4), atom_pairs in bonded[mol]['BON'].get('QUA', {}).items():
            for at1_id, at2_id in atom_pairs:
                q1, q2 = _bond_pair_qualified(mol, at1_id, at2_id, atom_map)
                if q1 is None:
                    continue
                key = (tuple(sorted([q1, q2])), (r0, k2, k3, k4))
                if key not in written:
                    written.add(key)
                    entries.append((q1, q2, r0, k2, k3, k4))

        for (r0, k2), atom_pairs in bonded[mol]['BON'].get('HAR', {}).items():
            for at1_id, at2_id in atom_pairs:
                q1, q2 = _bond_pair_qualified(mol, at1_id, at2_id, atom_map)
                if q1 is None:
                    continue
                key = (tuple(sorted([q1, q2])), (r0, k2, 0.0, 0.0))
                if key not in written:
                    written.add(key)
                    entries.append((q1, q2, r0, k2, 0.0, 0.0))

    if not entries:
        return ''

    lines = ['<CustomBondForce energy="(k2/2)*(r-r0)^2 + (k3/3)*(r-r0)^3 + (k4/4)*(r-r0)^4">',
             '<PerBondParameter name="r0"/>', '<PerBondParameter name="k2"/>',
             '<PerBondParameter name="k3"/>', '<PerBondParameter name="k4"/>']
    for q1, q2, r0, k2, k3, k4 in entries:
        r0_nm = r0 * 0.1
        k2_kj = k2 * 418.4
        k3_kj = k3 * 4184.0
        k4_kj = k4 * 41840.0
        lines.append(f'<Bond class1="{q1}" class2="{q2}" r0="{r0_nm}" k2="{k2_kj}" k3="{k3_kj}" k4="{k4_kj}"/>')
    lines.append('</CustomBondForce>')
    return '\n'.join(lines)


def _bond_pair_qualified(mol, at1_id, at2_id, atom_map):
    """Return (qual1, qual2) or (None, None) if either atom is absent/special."""
    a1 = atom_map.get(at1_id)
    a2 = atom_map.get(at2_id)
    if a1 is None or a2 is None:
        return None, None
    if a1[1] in ('NETF', 'TORQ') or a2[1] in ('NETF', 'TORQ'):
        return None, None
    return _qualify(mol, a1[1]), _qualify(mol, a2[1])


def gen_angle_force(bonded, mol_names):
    """``<CustomAngleForce>`` for HAR angles (theta0 in radians)."""
    entries, written = [], set()

    for mol in mol_names:
        if 'ANG' not in bonded[mol]:
            continue
        atom_map = _atom_map(bonded[mol]['ATO']['All'])

        for (theta0_deg, k), triples in bonded[mol]['ANG'].get('HAR', {}).items():
            theta0_rad = theta0_deg * math.pi / 180.0
            for at1_id, at2_id, at3_id in triples:
                a1 = atom_map.get(at1_id); a2 = atom_map.get(at2_id); a3 = atom_map.get(at3_id)
                if any(a is None or a[1] in ('NETF', 'TORQ') for a in (a1, a2, a3)):
                    continue
                q1 = _qualify(mol, a1[1]); q2 = _qualify(mol, a2[1]); q3 = _qualify(mol, a3[1])
                key = (tuple(sorted([q1, q3])) + (q2,), (theta0_deg, k))
                if key not in written:
                    written.add(key)
                    entries.append((q1, q2, q3, theta0_rad, k))

    if not entries:
        return ''

    lines = ['<CustomAngleForce energy="(k/2)*(theta-theta0)^2">',
             '<PerAngleParameter name="k"/>', '<PerAngleParameter name="theta0"/>']
    for q1, q2, q3, theta0, k in entries:
        k_kj = k * 4.184
        lines.append(f'<Angle class1="{q1}" class2="{q2}" class3="{q3}" theta0="{theta0}" k="{k_kj}"/>')
    lines.append('</CustomAngleForce>')
    return '\n'.join(lines)


def gen_dihedral_force(bonded, mol_names):
    """``<PeriodicTorsionForce>`` for NCO dihedrals.

    NCO tuple key in bonded dict: ``(k_kcal, periodicity, phase_rad)``;
    ``k`` is converted kcal/mol → kJ/mol, phase is already in radians.
    """
    entries, written = [], set()

    for mol in mol_names:
        if 'DIH' not in bonded[mol]:
            continue
        atom_map = _atom_map(bonded[mol]['ATO']['All'])

        for (k_kcal, periodicity, phase_rad), quads in bonded[mol]['DIH'].get('NCO', {}).items():
            k_kj = k_kcal * 4.184
            for at1_id, at2_id, at3_id, at4_id in quads:
                atoms = [atom_map.get(aid) for aid in (at1_id, at2_id, at3_id, at4_id)]
                if any(a is None or a[1] in ('NETF', 'TORQ') for a in atoms):
                    continue
                q = [_qualify(mol, a[1]) for a in atoms]
                key  = tuple(q) + (k_kcal, periodicity, phase_rad)
                rkey = tuple(reversed(q)) + (k_kcal, periodicity, phase_rad)
                if key not in written and rkey not in written:
                    written.add(key)
                    entries.append((*q, int(periodicity), phase_rad, k_kj))

    if not entries:
        return ''

    lines = ['<PeriodicTorsionForce>']
    for q1, q2, q3, q4, period, phase, k in entries:
        lines.append(f'<Proper class1="{q1}" class2="{q2}" class3="{q3}" class4="{q4}" '
                     f'periodicity1="{period}" phase1="{phase}" k1="{k}"/>')
    lines.append('</PeriodicTorsionForce>')
    return '\n'.join(lines)


# --------------------------------------------------------------------------
# CustomNonbondedForce builders (matrix lookup tables)
# --------------------------------------------------------------------------

def gen_exp_force(exp_entries, atom_types, bond_cutoff=2):
    """``<CustomNonbondedForce>`` for EXP: ``U(r) = A*exp(-alpha*r)``.

    A: kcal/mol → kJ/mol (×4.184); alpha: Å⁻¹ → nm⁻¹ (×10).
    """
    qualified = [q for _, _, q in atom_types]
    n   = len(qualified)
    idx = {q: i for i, q in enumerate(qualified)}
    a_mat     = _matrix(n)
    alpha_mat = _matrix(n)

    for q1, q2, A, alpha in exp_entries:
        i, j = idx[q1], idx[q2]
        a_mat[i][j]     = a_mat[j][i]     = A * 4.184
        alpha_mat[i][j] = alpha_mat[j][i] = alpha * 10.0

    return _custom_nb_xml(
        energy='a*exp(-alpha*r); a=aTable(t1,t2); alpha=alphaTable(t1,t2)',
        tables=[('aTable', a_mat), ('alphaTable', alpha_mat)],
        qualified=qualified, n=n, bond_cutoff=bond_cutoff,
    )


def gen_srd_force(entries, power, atom_types, bond_cutoff=2):
    """``<CustomNonbondedForce>`` for SRD/POW: ``U(r) = disp/(r^|p| + r0^|p|)``.

    disp: kcal·Å^|p| → kJ·nm^|p|; r0: Å → nm.
    """
    qualified = [q for _, _, q in atom_types]
    n   = len(qualified)
    idx = {q: i for i, q in enumerate(qualified)}
    disp_mat = _matrix(n)
    r0_mat   = _matrix(n)

    for q1, q2, P1, r0 in entries:
        i, j = idx[q1], idx[q2]
        disp_mat[i][j] = disp_mat[j][i] = P1 * 4.184 * (0.1 ** (-power))
        r0_mat[i][j]   = r0_mat[j][i]   = r0 * 0.1

    abs_power = abs(power)
    ps = str(int(abs_power)) if abs_power == int(abs_power) else str(abs_power)
    energy = (f'disp/(r^{ps} + r0^{ps}); '
              'disp=dispTable(t1,t2); r0=rTable(t1,t2)')

    return _custom_nb_xml(
        energy=energy,
        tables=[('dispTable', disp_mat), ('rTable', r0_mat)],
        qualified=qualified, n=n, bond_cutoff=bond_cutoff,
    )


def gen_str_force(str_entries, atom_types, bond_cutoff=2):
    """``<CustomNonbondedForce>`` for STR/STRC (shifted-truncated power)."""
    qualified = [q for _, _, q in atom_types]
    n   = len(qualified)
    idx = {q: i for i, q in enumerate(qualified)}
    p1_mat = _matrix(n)
    p2_mat = _matrix(n)
    p3_mat = _matrix(n)

    for q1, q2, P1, P2, P3 in str_entries:
        i, j = idx[q1], idx[q2]
        p1_mat[i][j] = p1_mat[j][i] = P1 * 4.184 * (0.1 ** P2)
        p2_mat[i][j] = p2_mat[j][i] = P2
        p3_mat[i][j] = p3_mat[j][i] = P3 * 0.1

    # Avoid 0^0 in the energy expression for unused pairs.
    for i in range(n):
        for j in range(n):
            if p1_mat[i][j] == 0.0:
                if p2_mat[i][j] == 0.0:
                    p2_mat[i][j] = 1.0
                if p3_mat[i][j] == 0.0:
                    p3_mat[i][j] = 1.0

    energy = ('step(p3 - r) * p1 * ((1/r^p2) - (1/p3^p2) + (p2*(r - p3))/(p3^(p2 + 1))); '
              'p1=p1Table(t1,t2); p2=p2Table(t1,t2); p3=p3Table(t1,t2)')

    return _custom_nb_xml(
        energy=energy,
        tables=[('p1Table', p1_mat), ('p2Table', p2_mat), ('p3Table', p3_mat)],
        qualified=qualified, n=n, bond_cutoff=bond_cutoff,
    )


def _custom_nb_xml(energy, tables, qualified, n, bond_cutoff=2):
    """Build a ``<CustomNonbondedForce>`` block with N×N Discrete2D lookup tables."""
    lines = [f'<CustomNonbondedForce energy="{energy}" bondCutoff="{bond_cutoff}">']
    for name, mat in tables:
        lines.append(f'<Function name="{name}" type="Discrete2D" xsize="{n}" ysize="{n}">')
        lines.append(_matrix_str(mat))
        lines.append('</Function>')
    lines.append('<PerParticleParameter name="t"/>')
    for i, q in enumerate(qualified):
        lines.append(f'<Atom class="{q}" t="{i}"/>')
    lines.append('</CustomNonbondedForce>')
    return '\n'.join(lines)


# --------------------------------------------------------------------------
# Charge penetration
# --------------------------------------------------------------------------

#: Coulomb constant in the .off's units, kcal·Å/mol/e². Same value pycryoff's CPN term uses
#: (``openmm_backend.term_builder``), so the two evaluate the identical function.
FCOV_COU = 332.0637


def gen_cpn_force(cpn_entries, atom_types, bond_cutoff=2):
    """``<CustomNonbondedForce>`` for CPN, charge penetration on the core/valence split.

    pycryoff writes the correction as (r in Å, energy in kcal/mol)::

        U = -FcovCou*scale/r * ( ZiQj*(1 + bj*r/2)*exp(-bj*r)
                               + ZjQi*(1 + bi*r/2)*exp(-bi*r)
                               + (V0 + V1*r + V2*r^2 + V3*r^3)*exp(-bi*r)
                               + (W0 + W1*r)*exp(-bj*r) )

    **This form is not symmetric in i and j, and that is the whole difficulty.** A
    ``CustomNonbondedForce`` decides for itself which particle of a pair is ``t1``; if the
    lookup tables are filled as though ``(i, j)`` and ``(j, i)`` were the same entry, half
    the pairs in the system silently evaluate the correction with the two atoms' valence
    widths swapped. It is a plausible-looking, wrong number, on exactly the pairs -- unlike
    types -- that carry the largest correction.

    So the expression is rewritten as a sum over the two *exponentials* rather than over the
    two atoms, with the ``(1 + b*r/2)`` factor folded into each polynomial (it always shares
    its own side's ``b``)::

        U = -( exp(-bA*r)*(a0 + a1*r + a2*r^2 + a3*r^3)
             + exp(-bB*r)*(c0 + c1*r + c2*r^2 + c3*r^3) ) / r

        side A (exponent bi):  a = [ ZjQi + V0, ZjQi*bi/2 + V1, V2, V3 ]
        side B (exponent bj):  c = [ ZiQj + W0, ZiQj*bj/2 + W1,  0,  0 ]

    Now exchanging the two atoms exchanges side A with side B and nothing else, so entry
    ``(j, i)`` is entry ``(i, j)`` with the two sides swapped -- and the total is the same
    whichever way OpenMM orders the pair. The B-side polynomial is padded to a cubic purely
    so that the swap is an exchange of identical shapes.

    Units: the ten table values absorb ``FcovCou``, the fitted ``scale``, kcal → kJ (×4.184)
    and Å → nm, leaving ``b`` in nm⁻¹ and the energy in kJ/mol.

    Returns ``''`` when there are no CPN pairs.
    """
    if not cpn_entries:
        return ''

    qualified = [q for _, _, q in atom_types]
    n   = len(qualified)
    idx = {q: i for i, q in enumerate(qualified)}
    names = ('bA', 'a0', 'a1', 'a2', 'a3', 'bB', 'c0', 'c1', 'c2', 'c3')
    mats  = {name: _matrix(n) for name in names}

    for q1, q2, params in cpn_entries:
        scale, ziqj, zjqi, bi, bj, v0, v1, v2, v3, w0, w1 = (float(v) for v in params)
        pre = 4.184 * FCOV_COU * scale
        # (coefficient of r^k) * 10^(k-1): the 1/r outside and r^k inside both change units.
        side_a = [pre * q * 10.0 ** (k - 1)
                  for k, q in enumerate((zjqi + v0, zjqi * bi / 2.0 + v1, v2, v3))]
        side_b = [pre * q * 10.0 ** (k - 1)
                  for k, q in enumerate((ziqj + w0, ziqj * bj / 2.0 + w1, 0.0, 0.0))]
        forward = (bi * 10.0, *side_a, bj * 10.0, *side_b)
        swapped = (bj * 10.0, *side_b, bi * 10.0, *side_a)
        i, j = idx[q1], idx[q2]
        # Swapped first, so a like pair (i == j, where the two orientations are the same
        # interaction anyway) ends up holding the deck's own orientation rather than its
        # mirror. Only tidiness: the energy is identical either way when bi == bj.
        for name, fwd, swp in zip(names, forward, swapped):
            mats[name][j][i] = swp
            mats[name][i][j] = fwd

    energy = ('-(exp(-bA*r)*(a0 + r*(a1 + r*(a2 + r*a3)))'
              ' + exp(-bB*r)*(c0 + r*(c1 + r*(c2 + r*c3))))/r; '
              + '; '.join(f'{name}={name}Table(t1,t2)' for name in names))

    return _custom_nb_xml(
        energy=energy,
        tables=[(f'{name}Table', mats[name]) for name in names],
        qualified=qualified, n=n, bond_cutoff=bond_cutoff,
    )


# --------------------------------------------------------------------------
# Explicit mutual polarization
# --------------------------------------------------------------------------

#: Å³ → nm³. Polarizabilities are quoted in Å³ everywhere in the .off/JSON world.
ANGSTROM3_TO_NM3 = 1.0e-3


def gen_multipole_force(bonded, mol_names, atom_types, type_names, type_to_charge,
                        polarization):
    """``<AmoebaMultipoleForce>`` carrying the permanent charges and the induced dipoles.

    A ``[POL]`` deck's electrostatics is *not* a set of point charges. It is those charges
    plus the dipoles they mutually induce, and induction is the only many-body term in the
    model -- every other form is a function of one pair distance. Writing such a fit as a
    plain ``<NonbondedForce>`` produces a different force field that looks like the right
    one, which is why :mod:`afmtogmx.core.read_json` warns when it sees the card.

    Dipoles and quadrupoles are all zero: this is a point-charge model *with* polarizability,
    not a multipole one. ``AmoebaMultipoleForce`` is used because it is OpenMM's only mutual
    induction solver, and with ``kz = kx = 0`` every site takes ``NoAxisType`` and needs no
    local frame -- which is what makes the zero higher multipoles well defined.

    Two conventions have to survive the trip, and only one of them fits in the XML:

    * **Polarization groups.** ``pgrp`` links each type to the types it is bonded to, so the
      group closes over the molecule and a molecule does not polarize itself. That is the
      standard AMOEBA choice and the one distributed atomic polarizabilities are derived
      under; without it an isolated molecule acquires induced dipoles out of nothing.
    * **Permanent intramolecular electrostatics.** ``INTRA=EXCLUDE`` switches off *all* of
      it. OpenMM's covalent maps come from the topology and scale 1-2 and 1-3 by zero but
      1-4 by 0.4 and 1-5 by 0.8, so any molecule wider than three bonds keeps a fraction of
      a term the fit did not have. There is no XML attribute for this -- the scale factors
      are not settable -- so it must be corrected on the ``System`` afterwards. That is what
      ``afm_openmm.match_multipole_covalent_maps`` does, and ``prepare_afm_system`` refuses
      to leave an ``AmoebaMultipoleForce`` untouched.

    ``poltype`` and the convergence tolerance are ``createSystem`` keyword arguments in
    OpenMM, not attributes of the XML element, so a deck that is not the OpenMM defaults
    (``mutual``, 1e-5) raises rather than being written out as something else.
    """
    if not polarization:
        return ''
    alphas = polarization.get('alphas') or {}
    thole  = float(polarization.get('thole', 0.39))

    poltype = str(polarization.get('poltype', 'mutual')).lower()
    if poltype != 'mutual':
        raise ValueError(
            f"[POL] POLTYPE={poltype.upper()} cannot be expressed in a ForceField XML: "
            "OpenMM takes the polarization type as a createSystem(polarization=...) "
            "argument, and the file would silently be read as MUTUAL.")
    eps = polarization.get('eps')
    if eps is not None and abs(float(eps) - 1e-5) > 1e-12:
        raise ValueError(
            f"[POL] EPS={eps} is not OpenMM's default 1e-5 and the XML has no attribute "
            "for it; pass createSystem(mutualInducedTargetEpsilon=...) and remove this "
            "check only once the caller does.")

    names = type_names or build_type_names(atom_types)
    raw_of = {qualified: raw for _, raw, qualified in atom_types}

    # pgrp entries are declared per TYPE, over bonds, so build type -> bonded types.
    neighbours = defaultdict(set)
    for mol in mol_names:
        atom_map = _atom_map(bonded[mol]['ATO']['All'])
        for pairs_dict in bonded[mol].get('BON', {}).values():
            for atom_pairs in pairs_dict.values():
                for at1_id, at2_id in atom_pairs:
                    q1, q2 = _bond_pair_qualified(mol, at1_id, at2_id, atom_map)
                    if q1 is None:
                        continue
                    neighbours[q1].add(q2)
                    neighbours[q2].add(q1)

    lines = ['<AmoebaMultipoleForce>']
    for _, raw, qualified in atom_types:
        c0 = type_to_charge.get(qualified, 0.0)
        lines.append(f'<Multipole type="{names[qualified]}" kz="0" kx="0" c0="{c0}" '
                     'd1="0.0" d2="0.0" d3="0.0" '
                     'q11="0.0" q21="0.0" q22="0.0" q31="0.0" q32="0.0" q33="0.0"/>')
    for _, raw, qualified in atom_types:
        alpha = float(alphas.get(raw, alphas.get(qualified, 0.0))) * ANGSTROM3_TO_NM3
        group = sorted((names[other] for other in neighbours.get(qualified, ())
                        if other in names), key=_type_sort_key)
        pgrp = ''.join(f' pgrp{k}="{name}"' for k, name in enumerate(group, start=1))
        lines.append(f'<Polarize type="{names[qualified]}" polarizability="{alpha}" '
                     f'thole="{thole}"{pgrp}/>')
    lines.append('</AmoebaMultipoleForce>')
    return '\n'.join(lines)


# --------------------------------------------------------------------------
# Exclusions
# --------------------------------------------------------------------------

def required_bond_cutoff(bonded, mol_names, max_cutoff=5):
    """The ``bondCutoff`` that reproduces the deck's ``[Exc]`` card exactly.

    OpenMM cannot be handed an exclusion list for a ``CustomNonbondedForce`` in a ForceField
    XML; it takes a ``bondCutoff`` and derives the pairs from the topology. The deck states
    its exclusions explicitly instead, and the two only agree by coincidence. Hard-coding 2
    is right for a deck that excludes 1-2 and 1-3, and wrong for the hybrid decks used here,
    which exclude every intramolecular pair -- leaving acetonitrile's three 1-4 H···N pairs
    interacting through an exchange wall fitted to contact distances, and through the
    penetration correction, neither of which the fit contained.

    So derive it, and refuse to guess: if no cutoff reproduces the card, raise.
    """
    needed = 0
    for mol in mol_names:
        atom_map = _atom_map(bonded[mol]['ATO']['All'])
        real = {aid for aid, (_, attype) in atom_map.items()
                if attype not in ('NETF', 'TORQ')}

        adj = defaultdict(set)
        for pairs_dict in bonded[mol].get('BON', {}).values():
            for atom_pairs in pairs_dict.values():
                for a, b in atom_pairs:
                    if a in real and b in real:
                        adj[a].add(b)
                        adj[b].add(a)

        # Bond-graph distance between every pair of real atoms.
        dist = {}
        for src in real:
            seen = {src: 0}
            frontier = [src]
            while frontier:
                nxt = []
                for a in frontier:
                    for b in adj[a]:
                        if b not in seen:
                            seen[b] = seen[a] + 1
                            nxt.append(b)
                frontier = nxt
            for b, d in seen.items():
                if b != src:
                    dist[frozenset((src, b))] = d

        declared = {frozenset((int(a), int(b))) for a, b in bonded[mol].get('EXC', [])
                    if int(a) in real and int(b) in real}

        for cutoff in range(0, max_cutoff + 1):
            if {pair for pair, d in dist.items() if d <= cutoff} == declared:
                needed = max(needed, cutoff)
                break
        else:
            raise ValueError(
                f"{mol}: the [Exc] card ({len(declared)} intramolecular pair(s)) is not the "
                f"set any bondCutoff up to {max_cutoff} produces from the molecule's bonds. "
                "An OpenMM CustomNonbondedForce in a ForceField XML can only exclude by bond "
                "count, so the exported model would not be the fitted one.")
    return needed


# --------------------------------------------------------------------------
# File output
# --------------------------------------------------------------------------

def write_xml(path, sections):
    """Wrap ``sections`` in ``<ForceField>`` and write to ``path``."""
    with open(path, 'w') as f:
        f.write('<ForceField>\n\n')
        f.write('\n\n'.join(sections) + '\n')
        f.write('\n</ForceField>\n')
