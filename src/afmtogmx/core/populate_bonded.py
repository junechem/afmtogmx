"""Build a `bonded[mol]` dict for a new molecule from a user-supplied topology
directory.

The directory holds four files:

  atoms.dat            atom_number  atom_type
  bonds.dat            atom_a       atom_b
  valid_dihedrals.dat  type1-type2-type3-type4   (one signature per line)
  parameters.dat       #define <name> <subtype> <values...>

This module is a parameter-binding layer: parsers produce primitive structures,
graph traversal derives bonds/angles/dihedrals from connectivity, and the
orchestrator returns a `bonded[mol]` dict shaped exactly like
:func:`afmtogmx.core.functions._parse_bonded` produces.  The orchestrator does
not mutate any `ReadOFF` state — `ReadOFF.populate_bonded` does that.
"""

import os

from afmtogmx.core import functions


# ---------------------------------------------------------------------------
# File parsers
# ---------------------------------------------------------------------------

def _parse_atoms(path):
    """Parse `atoms.dat`. Returns dict mapping atom number (int) -> atom type (str)."""
    atoms = {}
    with open(path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(
                    f"{path}: expected 'atom_number atom_type', got: {raw!r}"
                )
            atoms[int(parts[0])] = parts[1]
    return atoms


def _parse_bonds(path):
    """Parse `bonds.dat`. Returns list of (a, b) atom-number tuples."""
    bonds = []
    with open(path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(
                    f"{path}: expected 'atom_a atom_b', got: {raw!r}"
                )
            bonds.append((int(parts[0]), int(parts[1])))
    return bonds


def _parse_valid_dihedrals(path):
    """Parse `valid_dihedrals.dat`. Returns set of canonical type-signature strings."""
    valid = set()
    with open(path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            valid.add(_canonical_dihedral_sig(line))
    return valid


def _parse_parameters(path):
    """Parse `parameters.dat` `#define` lines.

    Returns a dict keyed by underscore count (1 = bond, 2 = angle, 3 = dihedral)
    whose value is `{signature: {'subtype': str, 'values': tuple[float, ...]}}`.
    """
    parameters = {1: {}, 2: {}, 3: {}}
    with open(path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line.startswith('#define'):
                continue
            parts = line.split()
            if len(parts) < 4:
                raise ValueError(
                    f"{path}: malformed #define (need name, subtype, values): {raw!r}"
                )
            name = parts[1]
            subtype = parts[2]
            try:
                values = tuple(float(v) for v in parts[3:])
            except ValueError as e:
                raise ValueError(
                    f"{path}: non-numeric parameter value in {raw!r}: {e}"
                )
            n_underscores = name.count('_')
            if n_underscores not in parameters:
                raise ValueError(
                    f"{path}: signature '{name}' has {n_underscores} underscores; "
                    f"expected 1 (bond), 2 (angle), or 3 (dihedral)"
                )
            parameters[n_underscores][name] = {'subtype': subtype, 'values': values}
    return parameters


# ---------------------------------------------------------------------------
# Graph traversal
# ---------------------------------------------------------------------------

def _build_adjacency(bonds):
    """Build undirected adjacency list from `(a, b)` bond tuples."""
    adj = {}
    for a, b in bonds:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    return adj


def _find_angles(adjacency):
    """Enumerate angles `(i, j, k)` where `j` is the central atom.

    Each unordered triple is emitted once with `i < k` (so `(i, j, k)` and
    `(k, j, i)` are not both returned).
    """
    angles = []
    for central, neighbors in adjacency.items():
        if len(neighbors) < 2:
            continue
        sorted_neighbors = sorted(neighbors)
        for i in range(len(sorted_neighbors)):
            for j in range(i + 1, len(sorted_neighbors)):
                angles.append((sorted_neighbors[i], central, sorted_neighbors[j]))
    return angles


def _find_dihedrals(angles, adjacency):
    """Enumerate unique dihedrals `(i, j, k, l)` from each angle by extending
    one bond on either side.

    A dihedral and its reverse are deduplicated.
    """
    seen = set()
    dihedrals = []

    def _add(quad):
        canonical = min(quad, quad[::-1])
        if canonical in seen:
            return
        seen.add(canonical)
        dihedrals.append(quad)

    for (a, b, c) in angles:
        for w in adjacency[a]:
            if w == b:
                continue
            _add((w, a, b, c))
        for x in adjacency[c]:
            if x == b:
                continue
            _add((a, b, c, x))
    return dihedrals


# ---------------------------------------------------------------------------
# Atom-type grouping (canonical signatures)
# ---------------------------------------------------------------------------

def _canonical_bond_sig(t1, t2):
    """Alphabetically sort the two atom types."""
    return '_'.join(sorted([t1, t2]))


def _canonical_angle_sig(t1, t2, t3):
    """Sort terminal atom types only; central is fixed."""
    if t1 <= t3:
        return f'{t1}_{t2}_{t3}'
    return f'{t3}_{t2}_{t1}'


def _canonical_dihedral_sig(sig):
    """Canonicalize a dash- or underscore-separated dihedral type signature.

    Returns the lexicographically smaller of the forward and reversed forms,
    using underscores as the internal separator.
    """
    parts = sig.replace('-', '_').split('_')
    if len(parts) != 4:
        raise ValueError(f"dihedral signature must have 4 atom types: {sig!r}")
    forward = '_'.join(parts)
    reverse = '_'.join(parts[::-1])
    return min(forward, reverse)


def _group_bonds_by_type(bonds, atoms_by_num):
    grouped = {}
    for (a, b) in bonds:
        sig = _canonical_bond_sig(atoms_by_num[a], atoms_by_num[b])
        grouped.setdefault(sig, []).append([a, b])
    return grouped


def _group_angles_by_type(angles, atoms_by_num):
    grouped = {}
    for (a, b, c) in angles:
        sig = _canonical_angle_sig(atoms_by_num[a], atoms_by_num[b], atoms_by_num[c])
        grouped.setdefault(sig, []).append([a, b, c])
    return grouped


def _group_dihedrals_by_type(dihedrals, atoms_by_num, valid_signatures):
    """Group dihedrals by canonical signature, dropping any whose signature is
    not in `valid_signatures`."""
    grouped = {}
    for (a, b, c, d) in dihedrals:
        sig = _canonical_dihedral_sig(
            f'{atoms_by_num[a]}_{atoms_by_num[b]}_{atoms_by_num[c]}_{atoms_by_num[d]}'
        )
        if sig not in valid_signatures:
            continue
        grouped.setdefault(sig, []).append([a, b, c, d])
    return grouped


# ---------------------------------------------------------------------------
# Parameter binding
# ---------------------------------------------------------------------------
#
# `parameters.dat` values are written in GROMACS .top output convention
# (length: nm, energy: kJ/mol, dihedral parameter order: phi K mult).
# The in-memory `bonded[mol]` dict, by contrast, stores the raw .off-file
# values (Å, kcal/mol, NCO order: K mult phi). Topology generation in
# topology.py applies the .off → GROMACS unit conversion when writing the .top
# file, so we must store .off-unit values here — otherwise the conversion
# would be applied on top of already-GROMACS values, producing wrong output.
#
# The converters below translate parameters.dat (subtype, values) tuples into
# the in-memory representation. Conversion factors are exact: 1 nm = 10 Å,
# 1 kcal/mol = 4.184 kJ/mol.

_KCAL_PER_KJ = 1.0 / 4.184


def _convert_bond_params(subtype, values):
    if subtype == 'HAR':
        # both same order: (r0, K). nm → Å (×10); kJ/mol/nm² → kcal/mol/Å² (÷418.4)
        return (values[0] * 10.0, values[1] * _KCAL_PER_KJ * 0.01)
    raise ValueError(f"populate_bonded: unsupported bond subtype {subtype!r}")


def _convert_angle_params(subtype, values):
    if subtype == 'HAR':
        # both same order: (theta0_deg, K). degrees unchanged; kJ → kcal (÷4.184)
        return (values[0], values[1] * _KCAL_PER_KJ)
    raise ValueError(f"populate_bonded: unsupported angle subtype {subtype!r}")


def _convert_dihedral_params(subtype, values):
    if subtype == 'NCO' or subtype == 'COS':
        # parameters.dat order: (phi_deg, K_kJ, mult)
        # in-memory order:     (K_kcal,  mult,  phi_deg)  — see functions._parse_bonded_section
        return (values[1] * _KCAL_PER_KJ, values[2], values[0])
    if subtype == 'HAR':
        # both same order: (phi_deg, K). degrees unchanged; kJ → kcal
        return (values[0], values[1] * _KCAL_PER_KJ)
    raise ValueError(f"populate_bonded: unsupported dihedral subtype {subtype!r}")


_KIND_CONVERTERS = {
    'BON': _convert_bond_params,
    'ANG': _convert_angle_params,
    'DIH': _convert_dihedral_params,
}


def _lookup_param(signature, param_table):
    """Look up `signature` (and its reversed form) in `param_table`.

    Returns the entry or `None` if not found.
    """
    if signature in param_table:
        return param_table[signature]
    reversed_sig = '_'.join(signature.split('_')[::-1])
    if reversed_sig in param_table:
        return param_table[reversed_sig]
    return None


def _build_bonded_section(grouped, param_table, kind, missing):
    """Bind parameters onto grouped atom-tuple lists, converting GROMACS-unit
    values from parameters.dat into .off-unit in-memory keys.

    Parameters
    ----------
    grouped : dict[str, list[list[int]]]
        Signature -> list of atom-number tuples (as built by `_group_*_by_type`).
    param_table : dict[str, dict]
        `parameters[n_underscores]` for the appropriate kind.
    kind : str
        ``'BON'``, ``'ANG'``, or ``'DIH'`` — selects the unit converter.
    missing : list[str]
        Mutated in-place with any signatures that fail to resolve.

    Returns
    -------
    dict[str, dict[tuple, list[list[int]]]]
        Maps subtype string (e.g. ``'HAR'``, ``'NCO'``) to a parameter-keyed
        dict of atom-tuple lists, matching the existing data-model convention.
    """
    convert = _KIND_CONVERTERS[kind]
    section = {}
    for signature, atom_tuples in grouped.items():
        entry = _lookup_param(signature, param_table)
        if entry is None:
            missing.append(signature)
            continue
        subtype = entry['subtype']
        param_key = convert(subtype, entry['values'])
        subtype_dict = section.setdefault(subtype, {})
        subtype_dict.setdefault(param_key, []).extend(atom_tuples)
    return section


# ---------------------------------------------------------------------------
# Exclusions
# ---------------------------------------------------------------------------

def _generate_exclusions(adjacency, atom_numbers):
    """Generate 1-2 and 1-3 exclusions from the bond graph.

    For each atom `i` in sorted `atom_numbers`, emit one line listing every
    `j > i` reachable in 1 or 2 bonds. Atoms with no exclusions are omitted.
    """
    exclusions = []
    for i in sorted(atom_numbers):
        neighbors = set(adjacency.get(i, []))
        next_neighbors = set()
        for n in neighbors:
            for nn in adjacency.get(n, []):
                if nn != i:
                    next_neighbors.add(nn)
        excluded = sorted(j for j in (neighbors | next_neighbors) if j > i)
        if excluded:
            exclusions.append([i, *excluded])
    return exclusions


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def build_new_molecule_bonded(directory, atom_type_universe):
    """Read all four files in `directory` and return a fully populated
    `bonded[mol]` dict (ATO/BON/ANG/BD3/DIH/CDI/EXC).

    Parameters
    ----------
    directory : str
        Path to a directory containing ``atoms.dat``, ``bonds.dat``,
        ``valid_dihedrals.dat``, and ``parameters.dat``.
    atom_type_universe : set[str]
        Atom types known to the surrounding force field. Every type in
        ``atoms.dat`` must be a member.

    Returns
    -------
    dict
        Single-molecule bonded dict, ready to be assigned to
        ``ReadOFF.bonded[new_mol_name]``.

    Raises
    ------
    ValueError
        If any required file is missing, if an atom type in ``atoms.dat`` is
        not in ``atom_type_universe``, or if any bond/angle/dihedral signature
        cannot be resolved against ``parameters.dat``. All gaps are reported
        in a single error message.
    """
    paths = {
        'atoms': os.path.join(directory, 'atoms.dat'),
        'bonds': os.path.join(directory, 'bonds.dat'),
        'valid_dihedrals': os.path.join(directory, 'valid_dihedrals.dat'),
        'parameters': os.path.join(directory, 'parameters.dat'),
    }
    missing_files = [p for p in paths.values() if not os.path.isfile(p)]
    if missing_files:
        raise ValueError(
            "populate_bonded: missing required file(s):\n  "
            + "\n  ".join(missing_files)
        )

    atoms_by_num = _parse_atoms(paths['atoms'])
    bonds = _parse_bonds(paths['bonds'])
    valid_dihedral_sigs = _parse_valid_dihedrals(paths['valid_dihedrals'])
    parameters = _parse_parameters(paths['parameters'])

    unknown_types = sorted(
        set(atoms_by_num.values()) - set(atom_type_universe)
    )

    adjacency = _build_adjacency(bonds)
    angles = _find_angles(adjacency)
    dihedrals = _find_dihedrals(angles, adjacency)

    grouped_bonds = _group_bonds_by_type(bonds, atoms_by_num)
    grouped_angles = _group_angles_by_type(angles, atoms_by_num)
    grouped_dihedrals = _group_dihedrals_by_type(
        dihedrals, atoms_by_num, valid_dihedral_sigs
    )

    missing_bonds = []
    missing_angles = []
    missing_dihedrals = []

    BON = _build_bonded_section(grouped_bonds, parameters[1], 'BON', missing_bonds)
    ANG = _build_bonded_section(grouped_angles, parameters[2], 'ANG', missing_angles)
    DIH = _build_bonded_section(grouped_dihedrals, parameters[3], 'DIH', missing_dihedrals)

    # Every signature listed in valid_dihedrals.dat must resolve to a #define,
    # even if no dihedral actually appears for it in the graph.
    for sig in sorted(valid_dihedral_sigs):
        if _lookup_param(sig, parameters[3]) is None:
            if sig not in missing_dihedrals:
                missing_dihedrals.append(sig)

    if unknown_types or missing_bonds or missing_angles or missing_dihedrals:
        sections = []
        if unknown_types:
            sections.append(
                "  unknown atom types in atoms.dat (not present in any "
                "surviving molecule of the source ReadOFF):\n    "
                + ", ".join(unknown_types)
            )
        if missing_bonds:
            sections.append(
                "  bond signatures with no #define in parameters.dat:\n    "
                + ", ".join(sorted(set(missing_bonds)))
            )
        if missing_angles:
            sections.append(
                "  angle signatures with no #define in parameters.dat:\n    "
                + ", ".join(sorted(set(missing_angles)))
            )
        if missing_dihedrals:
            sections.append(
                "  dihedral signatures with no #define in parameters.dat:\n    "
                + ", ".join(sorted(set(missing_dihedrals)))
            )
        raise ValueError(
            "populate_bonded: cannot build bonded section. Fix all of the "
            "following and rerun:\n" + "\n".join(sections)
        )

    # Ensure every standard subtype slot exists, matching gen_empty_bonded(),
    # so downstream code that indexes BON['HAR'] etc. never KeyErrors.
    bonded = functions.gen_empty_bonded()
    bonded['ATO']['All'] = {n: (t, t) for n, t in atoms_by_num.items()}
    bonded['ATO']['Virtual'] = {}
    for subtype, params in BON.items():
        bonded['BON'].setdefault(subtype, {}).update(params)
    for subtype, params in ANG.items():
        bonded['ANG'].setdefault(subtype, {}).update(params)
    for subtype, params in DIH.items():
        bonded['DIH'].setdefault(subtype, {}).update(params)
    bonded['EXC'] = _generate_exclusions(adjacency, atoms_by_num.keys())

    return bonded
