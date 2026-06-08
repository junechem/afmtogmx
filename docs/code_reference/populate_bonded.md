# `populate_bonded.py` — build a molecule's bonded dict from a topology directory

Reconstructs a `bonded[mol]` dict for a new molecule from four user-supplied
`.dat` files, shaped exactly like `functions._parse_bonded` produces. Pure
parameter-binding layer: parsers → graph traversal → orchestrator. Does **not**
mutate any `ReadOFF` state — `ReadOFF.populate_bonded` (in `gen_md.py`) does that.

**Input directory files:**
```
atoms.dat            atom_number  atom_type
bonds.dat            atom_a       atom_b
valid_dihedrals.dat  type1-type2-type3-type4   (one signature per line)
parameters.dat       #define <name> <subtype> <values...>
```

**Unit convention:** `parameters.dat` values are in GROMACS `.top` output units
(nm, kJ/mol, dihedral order `phi K mult`). The in-memory `bonded` dict stores raw
`.off` units (Å, kcal/mol, NCO order `K mult phi`), because `topology.py` applies
the `.off`→GROMACS conversion on write. The `_convert_*` functions here translate
`parameters.dat` → in-memory `.off` units. Constant: `_KCAL_PER_KJ = 1/4.184`.

## File parsers
- `_parse_atoms(path)` — `{atom_number(int): atom_type(str)}`.
- `_parse_bonds(path)` — `[(a, b), …]` atom-number tuples.
- `_parse_valid_dihedrals(path)` — set of canonical dihedral signature strings.
- `_parse_parameters(path)` — `{underscore_count: {signature: {'subtype', 'values'}}}`
  where 1=bond, 2=angle, 3=dihedral. Validates `#define` lines.

## Graph traversal
- `_build_adjacency(bonds)` — undirected adjacency list `{atom: [neighbors]}`.
- `_find_angles(adjacency)` — angles `(i, j, k)` with `j` central, each unordered
  triple once (`i < k`).
- `_find_dihedrals(angles, adjacency)` — unique dihedrals `(i, j, k, l)` extending
  each angle by one bond; reverse-deduplicated.

## Canonical type signatures + grouping
- `_canonical_bond_sig(t1, t2)` — sorted `t1_t2`.
- `_canonical_angle_sig(t1, t2, t3)` — terminals sorted, central fixed.
- `_canonical_dihedral_sig(sig)` — smaller of forward/reverse, underscore-joined;
  accepts dash or underscore input. Raises if not 4 types.
- `_group_bonds_by_type` / `_group_angles_by_type` — `{signature: [atom_tuples]}`.
- `_group_dihedrals_by_type(dihedrals, atoms_by_num, valid_signatures)` — groups,
  dropping any signature not in `valid_signatures`.

## Parameter binding
- `_convert_bond_params(subtype, values)` — `HAR` only: nm→Å, kJ/mol/nm²→kcal/mol/Å².
- `_convert_angle_params(subtype, values)` — `HAR` only: degrees unchanged, kJ→kcal.
- `_convert_dihedral_params(subtype, values)` — `NCO`/`COS` reorder to `(K_kcal,
  mult, phi)`; `HAR` keeps order, kJ→kcal.
- `_KIND_CONVERTERS` — maps `'BON'/'ANG'/'DIH'` to the converters above.
- `_lookup_param(signature, param_table)` — looks up the signature and its reverse;
  returns the entry or `None`.
- `_build_bonded_section(grouped, param_table, kind, missing)` — binds parameters
  onto grouped atom tuples, converting units. Appends unresolved signatures to
  `missing`. Returns `{subtype: {param_key: [atom_tuples]}}`.

## Exclusions
- `_generate_exclusions(adjacency, atom_numbers)` — 1-2 and 1-3 exclusions; one line
  per atom `i` listing reachable `j > i`. Atoms with none are omitted.

## Orchestrator
- `build_new_molecule_bonded(directory, atom_type_universe)` — reads all four files
  and returns a fully populated single-molecule bonded dict (ATO/BON/ANG/BD3/DIH/
  CDI/EXC) starting from `functions.gen_empty_bonded()`. **Fails fast**, collecting
  *all* gaps into one `ValueError`: missing files, atom types in `atoms.dat` not in
  `atom_type_universe`, and any bond/angle/dihedral signature with no `#define`
  (including signatures listed in `valid_dihedrals.dat` that never resolve).
