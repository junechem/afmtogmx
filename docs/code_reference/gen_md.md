# `gen_md.py` — the `ReadOFF` class

Main entry point of the package. Parsing a `.off` file produces a `ReadOFF`
object that holds all force-field data and exposes two output backends
(`off.gmx`, `off.openmm`).

## Class `ReadOFF`

### Constructor
`__init__(self, off_loc)`
- Reads/parses the `.off` file at `off_loc` and populates instance attributes.
- Build sequence: `_gen_sections_dict()` → `_gen_bonded()` → build `self.charges`
  (all 0.0, keyed by atom name, NETF/TORQ excluded) → `_gen_nonbonded()` →
  build `self.residues` → instantiate `GROMACSBackend(self)` and
  `OpenMMBackend(self)`.

### Key attributes
- `off_loc` — path to the `.off` file.
- `_ff_bonded` — scratch dict for top-of-file bonding info (used during parse).
- `bonded` — `{molname: {ATO, BON, ANG, BD3, DIH, CDI, EXC}}`; fitted bonded params.
- `nonbonded` — `{(At1, At2): {InteractionType: [param_sets]}}`; atom pairs sorted.
- `charges` — `{molname: {atom_name: charge}}`; defaults 0.0.
- `residues` — `{'Definitions': {...}, 'Residues': {...}}` per molecule.
- `sections` — `{ff_input, intra_potential, inter_potential, molecular_definition, table_potential}` raw strings.
- `gmx` — `GROMACSBackend`; owns GROMACS config + tabpot outputs.
- `openmm` — `OpenMMBackend`; owns OpenMM config + XML/PDB outputs.

### Private parse methods
- `_gen_sections_dict()` — reads the file, splits into the 5 sections via
  `functions._find_off_keywords`, stores in `self.sections`. Raises on missing file.
- `_gen_bonded()` — recognizes keywords, filters into bonded/nonbonded, finds
  molnames, splits per molecule, gathers fitted bonded params, and calls
  `functions._parse_bonded` per molecule. Resets `functions.total_bonded_added` to 0.
- `_gen_nonbonded()` — cleans the inter-potential section and populates
  `self.nonbonded`; collapses any `COU*` interaction label to `'COU'`.

### Deprecated forwarders (emit `DeprecationWarning`)
All forward to the matching `off.gmx` member — kept for backward compatibility:
- Properties: `config`, `nonbonded_tabpot`, `bonded_tabpot` (getters + setters).
- Methods: `gen_nonbonded_tabpot`, `gen_bonded_tabpot`, `gen_nonbonded_topology`,
  `gen_bonded_topology`, `write_nonbonded_tabpot`, `write_bonded_tabpot`,
  `set_config`, `get_config`.
- **New code should call these on `off.gmx` directly.**

### Public methods (live on `ReadOFF`, not a backend)

`gen_residues(residue_definition={}, residue_atnums={})`
- Populates `self.residues` with custom residue groupings (by atom type and/or
  atom number). Validates against `self.bonded` via the `residues` module
  (`_check_residue_definitions`, `_check_residue_atnums`, `_set_*`). Prints progress.
- Formats: `residue_definition = {mol: {ResName: [AtType, ...]}}`;
  `residue_atnums = {mol: {ResName: [[atnum, ...], ...]}}`.

`change_molecule(mol_name, reference_ff, atom_name_map=None, ref_mol_name=None)`
- Replaces one molecule's bonded params + water-water nonbonded params with those
  from a stored reference FF in `src/afmtogmx/reference_ff/<reference_ff>.off`.
- Renames affected atom types via `atom_name_map`. **Cross-term pairs** (one atom
  in `mol_name`, one elsewhere) keep their original fitted values; only the
  `mol_name` side type name is renamed (those cross terms were fitted against the
  reference water model).
- Auto-loads a sibling `.charges` file if present. Rebuilds `self.residues`.
  Returns `self` for chaining.
- Raises `KeyError` (missing mol), `FileNotFoundError` (missing FF), `ValueError`
  (multi-molecule reference FF without `ref_mol_name`).

`populate_bonded(directory, new_mol_name, remove_molecules=None)`
- Returns a **new** `ReadOFF` (deepcopy; `self` unchanged) whose bonded section
  for `new_mol_name` is rebuilt from a topology directory via
  `populate_bonded.build_new_molecule_bonded`. Co-fitted molecules carry through;
  `remove_molecules` drops originals. Prunes `nonbonded` to surviving atom types.
  Adds default-0.0 charges and a residues entry for the new molecule.

`load_charges_from_file(file_path)`
- Reads a charge file into `self.charges`. Format: a lone molname line begins a
  block, then `atomname charge` lines. `#`/blank lines ignored. Unknown
  molecule/atom names warn and skip. An atom-charge pair seen before any molname
  is applied to every molecule containing that atom name. **Overwrites** existing
  charges for listed atoms. Returns `self`. Raises on missing file / parse error.
