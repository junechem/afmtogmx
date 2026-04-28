"""XML section builders for OpenMM force-field files.

Module-level helpers that turn ``ReadOFF`` data structures (``bonded``,
``nonbonded``, ``charges``) into OpenMM ``<ForceField>`` XML sections.
Used by :class:`afmtogmx.core.openmm_backend.OpenMMBackend`.

Atom types are *qualified* by molecule name ("UNK_C1") to keep them
globally unique even when two ``.off`` molecules happen to reuse the
same raw type label.
"""
import math
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
        element = attype[0]
        n = counters.get(element, 0)
        counters[element] = n + 1
        result[atid] = f'{element}{n}'
    return result


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
    """Return ``{qualified_type: charge}`` joining ATO atom names with charges dict."""
    name_lookup = {}  # (mol, attype) → atname for charge lookup
    for mol, raw, qualified in atom_types:
        for atnum, (atname, attype) in bonded[mol]['ATO']['All'].items():
            if attype == raw:
                name_lookup.setdefault((mol, raw), atname)

    result = {}
    for mol, raw, qualified in atom_types:
        atname = name_lookup.get((mol, raw))
        result[qualified] = charges.get(mol, {}).get(atname, 0.0)
    return result


def collect_nonbonded(nonbonded, atom_types):
    """Walk ``nonbonded`` and split entries by interaction type.

    Returns
    -------
    exp_entries  : list of ``(qual1, qual2, A, alpha)`` (kcal·Å units)
    str_entries  : list of ``(qual1, qual2, P1, P2, P3)`` (kcal·Å units)
    srd_by_power : dict ``{power: [(qual1, qual2, P1, r0)]}`` (kcal·Å units)

    Raw .off type names (``OW``) are expanded to all matching qualified
    types (``H2OQM_OW``, ``H2OMM_OW``); shared parameters get applied to
    every qualified pair, matching the .off file's intent.
    """
    raw_to_qualified = defaultdict(list)
    for mol, raw, qualified in atom_types:
        raw_to_qualified[raw].append(qualified)

    exp_entries, str_entries, srd_by_power = [], [], {}

    def _add_srd(q1, q2, P1, power, r0):
        srd_by_power.setdefault(power, []).append((q1, q2, P1, r0))

    for (raw1, raw2), params in nonbonded.items():
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
                        elif itype == 'BUC':
                            p1, p2, p3 = float(pset[0]), float(pset[1]), float(pset[2])
                            exp_entries.append((q1, q2, p1, p3))
                            _add_srd(q1, q2, p2, -6.0, 0.0)

    return exp_entries, str_entries, srd_by_power


# --------------------------------------------------------------------------
# Section builders
# --------------------------------------------------------------------------

def gen_atomtypes(bonded, atom_types):
    """``<AtomTypes>`` with one ``<Type>`` per qualified atom type."""
    lines = ['<AtomTypes>']
    for mol, raw, qualified in atom_types:
        if _is_virtual(bonded, mol, raw):
            lines.append(f'<Type name="{qualified}" class="{qualified}" mass="0.0"/>')
        else:
            element = raw[0]
            mass = _get_mass(element)
            lines.append(f'<Type name="{qualified}" class="{qualified}" element="{element}" mass="{mass}"/>')
    lines.append('</AtomTypes>')
    return '\n'.join(lines)


def gen_residues(bonded, mol_names, molname_translations):
    """``<Residues>`` listing atoms, bonds, and virtual sites for every molecule."""
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
            lines.append(f'<Atom name="{unique[atid]}" type="{_qualify(mol, attype)}"/>')

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

        lines.append('</Residue>')
    lines.append('</Residues>')
    return '\n'.join(lines)


def _virtual_site_xml(site_name, definition, unique):
    """Parse a virtual-site definition tuple → ``<VirtualSite>`` element."""
    try:
        parts     = list(definition)
        num_atoms = int(parts[0].rstrip(':'))
        weights, atom_ids = [], []
        i = 1
        while i < len(parts) and len(weights) < num_atoms:
            weights.append(float(parts[i]))
            atom_ids.append(int(parts[i + 1]))
            i += 2
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
    except (ValueError, IndexError):
        return None


def gen_nonbonded_force(atom_types, type_to_charge):
    """``<NonbondedForce>`` carrying point charges (sigma=epsilon=0)."""
    lines = ['<NonbondedForce coulomb14scale="1.0" lj14scale="1.0">']
    for mol, raw, qualified in atom_types:
        charge = type_to_charge.get(qualified, 0.0)
        lines.append(f'<Atom type="{qualified}" charge="{charge}" sigma="0.0" epsilon="0.0"/>')
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
        lines.append(f'<Bond type1="{q1}" type2="{q2}" r0="{r0}" k2="{k2}" k3="{k3}" k4="{k4}"/>')
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
        lines.append(f'<Angle type1="{q1}" type2="{q2}" type3="{q3}" theta0="{theta0}" k="{k}"/>')
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

def gen_exp_force(exp_entries, atom_types):
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
        qualified=qualified, n=n,
    )


def gen_srd_force(entries, power, atom_types):
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
        qualified=qualified, n=n,
    )


def gen_str_force(str_entries, atom_types):
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
        qualified=qualified, n=n,
    )


def _custom_nb_xml(energy, tables, qualified, n):
    """Build a ``<CustomNonbondedForce>`` block with N×N Discrete2D lookup tables."""
    lines = [f'<CustomNonbondedForce energy="{energy}" bondCutoff="2">']
    for name, mat in tables:
        lines.append(f'<Function name="{name}" type="Discrete2D" xsize="{n}" ysize="{n}">')
        lines.append(_matrix_str(mat))
        lines.append('</Function>')
    lines.append('<PerParticleParameter name="t"/>')
    for i, q in enumerate(qualified):
        lines.append(f'<Atom type="{q}" t="{i}"/>')
    lines.append('</CustomNonbondedForce>')
    return '\n'.join(lines)


# --------------------------------------------------------------------------
# File output
# --------------------------------------------------------------------------

def write_xml(path, sections):
    """Wrap ``sections`` in ``<ForceField>`` and write to ``path``."""
    with open(path, 'w') as f:
        f.write('<ForceField>\n\n')
        f.write('\n\n'.join(sections) + '\n')
        f.write('\n</ForceField>\n')
