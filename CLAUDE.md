# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Code Reference Docs (read these first)

`docs/code_reference/` contains a per-module reference describing every function in
`src/afmtogmx/core/` — its purpose, inputs/outputs, and non-obvious behavior. **Read
the relevant reference file before reading or modifying a `core/` source file**; it
is faster than re-reading the source and is intended to save that round trip. Start
with `docs/code_reference/README.md` for the module map and how the pieces fit
together. Covered modules: `gen_md`, `functions`, `tabulated_potentials`,
`topology`, `gmx_backend`, `openmm_backend`, `xml_generation`, `pdb_processing`,
`populate_bonded`, `report`.
(`compare.py`, `residues.py`, and `__init__.py` are not documented there.)

**Keep these docs in sync with the code.** Whenever you change, add, or remove a
function in a documented `core/` module, update the matching reference file in
`docs/code_reference/` in the same change. If a doc and the code disagree, the code
wins — fix the doc.

## Project Overview

`afmtogmx` is a Python package that converts .off files produced by CRYOFF into input files for GROMACS **and OpenMM** molecular dynamics simulations. For GROMACS it generates tabulated potentials and topology files for both bonded and nonbonded interactions; for OpenMM it generates a self-contained `<ForceField>` XML and preprocesses PDB files to match it. It can also rewrite the force field before output — swapping a molecule for a stored reference FF, or rebuilding a molecule's bonded terms from an atoms/bonds topology.

For the **authoritative meanings, parameter orders, and units of every CRYOFF `.off` interaction symbol** (HAR, QUA, BD3, DIH, CDIH, COU, THC, GLJ, BUC, DBU, STR, EXP, GEX, POW, CSP, GDP, PEX, DPO, SRD, TTP, etc.) and for known naming inconsistencies inside this package (`CDI` vs `CDIH`, `STR` vs `STRC` vs `SHTR`, `POW` vs `POW_6`), read `docs/CRYOFF_REFERENCE.md`. Anything in this CLAUDE.md that conflicts with that reference is wrong; fix it.

## Core Architecture

### Main Entry Point

The primary class is `ReadOFF` in `src/afmtogmx/core/gen_md.py`. All workflow operations start by instantiating this class:

```python
import afmtogmx as afm
off = afm.ReadOFF(off_loc="path/to/intra.off")
```

### Key Data Structures

After initialization, `ReadOFF` objects contain:
- `off.bonded`: Dictionary with fitted bonded parameters organized by molecule name and interaction type (ATO, BON, ANG, BD3, DIH, CDI, EXC)
- `off.nonbonded`: Dictionary with nonbonded parameters organized by atom pairs and interaction types
- `off.charges`: Dictionary with atomic charges per molecule (default 0.0 for all atoms; must be manually set if needed)
- `off.residues`: Dictionary with residue definitions and atom number mappings
- `off.sections`: Internal dictionary splitting .off file into 5 sections (ff_input, intra_potential, inter_potential, molecular_definition, table_potential)
- `off.gmx`: `GROMACSBackend` — GROMACS config plus `off.gmx.nonbonded_tabpot` / `off.gmx.bonded_tabpot` (None until the matching `gen_*` method runs)
- `off.openmm`: `OpenMMBackend` — OpenMM config and XML/PDB output methods

### Backend Architecture

Shared force-field data lives on `ReadOFF`; engine-specific output lives on a backend.

```
ReadOFF(off_loc)
  ├── off.bonded / off.nonbonded / off.charges / off.residues   ← shared parsed data
  ├── off.gmx      -> GROMACSBackend  (gmx_backend.py)
  │      ├── gen_nonbonded_tabpot / gen_bonded_tabpot   -> tabulated_potentials.py
  │      ├── write_nonbonded_tabpot / write_bonded_tabpot
  │      └── gen_nonbonded_topology / gen_bonded_topology -> topology.py
  └── off.openmm   -> OpenMMBackend   (openmm_backend.py)
         ├── gen_xml         -> xml_generation.py
         └── preprocess_pdb  -> pdb_processing.py
```

`ReadOFF.change_molecule`, `ReadOFF.populate_bonded` and `ReadOFF.gen_residues` rewrite
the shared data before either backend runs.

The old flat methods on `ReadOFF` (`off.set_config`, `off.gen_nonbonded_tabpot`, …) are
deprecated forwarders to `off.gmx` that emit `DeprecationWarning`. **Never write new code
against them, and never add examples using them to docs.**

### Configuration System

Each backend keeps its own config, so defaults are set once instead of passing them to
every method call. Both use the same parameter resolution order:

1. Explicitly passed parameter (if not None)
2. Configuration value
3. Built-in default

```python
off.gmx.set_config(spacing_nonbonded=0.001, scale_C6=False, tabpot_prefix='my_table')
off.gmx.get_config()               # full config dictionary
off.gmx.get_config('scale_C6')     # specific value
# set_config returns the backend, so calls chain
```

**`off.gmx` config keys** (`gmx_backend.py`):
- `incl_mol`: List of molecules to include (default: empty list, which means all)
- `excl_interactions`: List of interaction types to exclude (default: empty list)
- `excl_pairs`: List of atom pairs to exclude (default: empty list)
- `spacing_nonbonded`: Spacing for nonbonded tables (default: 0.0005 nm)
- `length_nonbonded`: Length for nonbonded tables (default: 3 nm)
- `scale_C6`: Enable C6 scaling for dispersion corrections (default: True)
- `sc_sigma`: Soft-core sigma for free energy calculations (default: 0.0)
- `spacing_bonded`: Spacing for bonded tables (default: 0.0001 nm)
- `length_bonded`: Length for bonded tables (default: 0.3 nm)
- `tabpot_prefix`: Prefix for tabulated potential files (default: 'table')
- `tabpot_dir`: Directory for tabulated potential files (default: '')
- `write_blank`: Write blank interactions (default: True)
- `name_translation`: **Atom** name translation dictionary (default: {})
- `special_pairs`: Custom attractive interactions per pair (default: empty dict)

**`off.openmm` config keys** (`openmm_backend.py`):
- `incl_mol`: List of molecules to include (default: empty list, which means all)
- `molname_translations`: **Molecule** name → PDB residue name (default: {}), e.g. `{'H2OQM': 'SOL'}`

**Three similarly-named things — do not confuse them**: `name_translation` (`off.gmx`
config, atom names), `molname_translations` (`off.openmm` config, plural, molecule → PDB
residue name), `molname_translation` (singular, a `write_report()` kwarg, display names).

### Module Responsibilities

**gen_md.py**
- Main `ReadOFF` class: parses .off files, populates the shared data structures,
  instantiates both backends
- Force-field-level public methods: `change_molecule`, `populate_bonded`,
  `gen_residues`, `load_charges_from_file`, `write_report`
- Deprecated forwarders to `off.gmx` (do not extend these)

**functions.py**
- Low-level parsing utilities for .off file format
- Interaction function definitions (exp, srd, shtr, powe, quarbond)
- Helper functions for filtering and organizing parameters
- Maintains global counter `total_bonded_added` for tracking

**gmx_backend.py**
- `GROMACSBackend`, exposed as `off.gmx`: GROMACS config plus the public GROMACS API
- Delegates to `tabulated_potentials.py` and `topology.py`

**openmm_backend.py**
- `OpenMMBackend`, exposed as `off.openmm`: OpenMM config plus `gen_xml` and `preprocess_pdb`
- Delegates to `xml_generation.py` and `pdb_processing.py`
- Imports OpenMM behind an availability guard — without it, XML masses fall back to 0.0

**tabulated_potentials.py**
- Generates tabulated potential tables for bonded and nonbonded interactions
- Handles special interaction types: POW_*, DPO_*, SRD_*, PEX_*, EXP, BUC, etc.
- C6 scaling for dispersion corrections
- Soft-core transformations for free energy calculations

**topology.py**
- Generates GROMACS .top topology files
- Writes [ nonbonded_params ] and [ pairs ] sections
- Template-based topology file modification
- Atom name translation support

**xml_generation.py**
- Builds the OpenMM `<ForceField>` XML sections from parsed FF data
- Namespaces atom types as `<MOLNAME>_<TYPE>`; converts kcal/mol→kJ/mol and Å→nm
- Supports a narrower nonbonded set than GROMACS: EXP, STR/STRC, SRD, POW (as SRD with
  r0=0), BUC (split into EXP + SRD)

**pdb_processing.py**
- Renames PDB atoms to match an XML's residue atom names and emits fresh CONECT records
- Reads only the XML, never in-memory state, so one XML can process many PDBs
- Topology (graph-isomorphism) matching when the PDB has CONECT records; positional
  matching with a `maxwarn` failsafe otherwise
- See `docs/PDB_PROCESSING.md`

**populate_bonded.py**
- Rebuilds a molecule's bonded dict from a user topology directory (atoms + bonds)
- Derives angles, dihedrals and 1-2/1-3 exclusions from the bond graph
- Fails fast, aggregating every gap into one `ValueError`

**report.py**
- `write_report`: fixed-width ASCII parameter summary in native .off kcal/Å units
- Column schemas live in the module-level `BONDED_SCHEMAS` dict

**residues.py**
- Manages residue definitions and atom number mappings
- Validates user-provided residue definitions

**compare.py**
- Compares two force fields (two ReadOFF objects)
- Generates difference strings for documentation

## Standard Workflows

### GROMACS

```python
import afmtogmx as afm

off = afm.ReadOFF(off_loc="path/to/intra.off")
off.load_charges_from_file('charges.txt')
off.gmx.set_config(tabpot_prefix='my_table', tabpot_dir='tables/', scale_C6=False)

# gen_* methods store their result on the backend; write_* methods pick it up
off.gmx.gen_nonbonded_tabpot()      # -> off.gmx.nonbonded_tabpot
off.gmx.gen_bonded_tabpot()         # -> off.gmx.bonded_tabpot
off.gmx.write_nonbonded_tabpot()
off.gmx.write_bonded_tabpot()

off.gmx.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')
off.gmx.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')
```

A custom dict can still be passed explicitly when generating variants:
```python
custom_tabpot = off.gmx.gen_nonbonded_tabpot(spacing_nonbonded=0.001)
off.gmx.write_nonbonded_tabpot(nonbonded_tabpot=custom_tabpot)
```

### OpenMM

```python
off = afm.ReadOFF(off_loc="path/to/intra.off")
off.load_charges_from_file('charges.txt')

off.openmm.set_config(molname_translations={'H2OQM': 'SOL'})
off.openmm.gen_xml(output='forcefield.xml')
off.openmm.preprocess_pdb('conf.pdb', 'forcefield.xml', 'conf_processed.pdb')
```

The XML is self-contained — no second force field file is merged at runtime.
`preprocess_pdb` reads only the XML, so one XML can process many PDBs.

See detailed workflow examples in `README.md` and `src/afmtogmx/__init__.py`.

## Force Field Modification

These methods rewrite the shared data and therefore affect both backends. Call them
before any `gen_*`.

### `change_molecule(mol_name, reference_ff, atom_name_map=None, ref_mol_name=None)`

Replaces one molecule's bonded params and its pure intra-molecule nonbonded pairs with
those from a stored reference FF in `src/afmtogmx/reference_ff/<reference_ff>.off`.

```python
off.change_molecule(mol_name='H2OQM', reference_ff='BLYPSP-4F', ref_mol_name='H2OQM')
```

- **Cross-term pairs keep their original fitted values**; only the `mol_name`-side type
  name is renamed. Those cross terms were fitted against the reference water model.
- Auto-loads a sibling `.charges` file if present; rebuilds `self.residues`; returns `self`.
- `ref_mol_name` is required when the reference FF defines more than one molecule
  (`BLYPSP-4F` does).
- `atom_name_map` (`{old_type: new_type}`) renames types at the same time.

### `populate_bonded(directory, new_mol_name, remove_molecules=None)`

Returns a **new** `ReadOFF` (deepcopy; `self` unchanged) whose bonded section for
`new_mol_name` is rebuilt from a topology directory of atoms and bonds. Angles,
dihedrals and 1-2/1-3 exclusions are derived from the bond graph.

Required files in `directory`:

| File | Format |
|------|--------|
| `atoms.dat` | `atom_number  atom_type` |
| `bonds.dat` | `atom_a  atom_b` |
| `valid_dihedrals.dat` | `type1-type2-type3-type4`, one per line (dash or underscore) |
| `parameters.dat` | `#define <signature> <subtype> <values...>` |

- Underscore count in a `parameters.dat` signature selects the term kind: 1=bond,
  2=angle, 3=dihedral.
- **`parameters.dat` values are in GROMACS units** (nm, kJ/mol, dihedral order
  `phi K mult`); the in-memory dict stores raw .off units (Å, kcal/mol, NCO order
  `K mult phi`). Conversion is automatic.
- `valid_dihedrals.dat` is a filter: unlisted dihedral signatures are dropped, and every
  listed signature must resolve to a `#define` even if no such dihedral exists.
- Supported subtypes: bonds `HAR`; angles `HAR`; dihedrals `NCO`, `COS`, `HAR`.
  `BD3` and `CDI` terms are **not** generated.
- Fails fast: every missing file, unknown atom type and unresolved signature is collected
  into one `ValueError`.
- New molecule charges start at 0.0; `nonbonded` is pruned to surviving atom types.

### `write_report(path, incl_mol=None, name_translation=None, molname_translation=None, notation='standard')`

Fixed-width ASCII summary of charges, bonded and nonbonded parameters, sized for Word at
default margins with Courier New 10pt. Units are the .off file's native kcal/Å — **no
conversion**. `notation` is `'standard'` (A, alpha, Cn, R0) or `'PN'` (the manual's
generic P1, P2, …); anything else raises `ValueError`. Nonbonded pairs are kept under a
molecule's section if *either* atom belongs to it, so cross terms survive `incl_mol`.

## Key Concepts

### Interaction Types

The package handles various interaction types with different naming conventions:
- **Variable power interactions**: POW, PEX, DPO, SRD (formatted as `TYPE_POWER`, e.g., `POW_6`)
- **Fixed interactions**: EXP, BUC, THC, COU (coulombic)
- **Default attractive interactions**: POW_6, DPO_6, SRD_6 (placed in columns 4-5 of tabulated potentials)
- **Repulsive interactions**: All others except COU (placed in columns 6-7)

### Special Pairs

Use `special_pairs` dictionary to customize which interactions are treated as attractive for specific atom pairs. This is incompatible with `scale_C6=True`.

Example: `special_pairs = {('At1', 'At2'): ['POW_6', 'POW_8']}`

### C6 Scaling

By default (`scale_C6=True`), columns 4-5 in tabulated potentials are scaled by the C6 parameter to enable GROMACS dispersion corrections. This assumes only one attractive interaction per pair.

### Manual Charge Assignment

Charges must be manually assigned to the `off.charges` dictionary. By default, all atom charges are 0.0. The dictionary format is:
```python
off.charges['MolName']['AtomName'] = charge_value
```

Example for water (TIP3P-like charges):
```python
off.charges['H20QM']['OQM'] = -0.82
off.charges['H20QM']['HQM'] = 0.41
off.charges['H20QM']['EQM'] = 0.0
```

**Alternative: Load charges from file**

Use the `load_charges_from_file()` method to read charges from a text file:

```python
off.load_charges_from_file('charges.txt')
```

File format:
```
MOLNAME1
Atom1 Charge1
Atom2 Charge2
MOLNAME2
Atom3 Charge3
```

- Blank lines and lines starting with `#` are ignored (comments)
- Any atoms not specified in the file remain at their default charge of 0.0
- WARNING: This method will overwrite any previously set charges for atoms specified in the file
- Returns `self` for method chaining
- Legacy format: an `atom charge` pair appearing **before any molecule name** is applied
  to every molecule containing that atom name

### Soft-Core Scaling

Use `sc_sigma` in `off.gmx.gen_nonbonded_tabpot()` for free energy calculations. Tables are scaled to work correctly with the sc-sigma value in GROMACS .mdp files.

### Recent Changes

**Backend architecture (latest)**:
- `ReadOFF` split into shared parsed data plus two backends, `off.gmx` (`GROMACSBackend`)
  and `off.openmm` (`OpenMMBackend`); each owns its own config
- The old flat methods on `ReadOFF` remain as forwarders to `off.gmx` that emit
  `DeprecationWarning` — do not write new code or docs against them

**OpenMM support**:
- `off.openmm.gen_xml()` writes a self-contained `<ForceField>` XML
- `off.openmm.preprocess_pdb()` renames PDB atoms to match the XML and emits CONECT
  records, via `pdb_processing.py`; topology matching is used when the PDB has CONECT
  records, positional matching otherwise

**Force field manipulation**:
- `ReadOFF.change_molecule()` replaces a molecule with a stored reference FF from
  `src/afmtogmx/reference_ff/` (ships `BLYPSP-4F`)
- `ReadOFF.populate_bonded()` rebuilds a molecule's bonded section from an atoms/bonds
  topology directory

**Reporting**:
- `ReadOFF.write_report()` writes a publication-style fixed-width parameter summary

**Earlier**:
- Configuration system with `set_config()` / `get_config()` and explicit → config →
  default parameter resolution
- `load_charges_from_file()` convenience method
- Tabulated potentials stored on the object, so `write_*` methods need no arguments
- `excl_pairs` filtering; C12 column scaling for HFE calculations
- All charge *calculation* functions removed — only manual assignment is supported

## Development Notes

### File Organization

- Core logic: `src/afmtogmx/core/`
- Package exports: `src/afmtogmx/__init__.py`
- Bundled reference force fields: `src/afmtogmx/reference_ff/` (`.off` plus optional `.charges`)
- Reference docs: `docs/` (`CRYOFF_REFERENCE.md`, `PDB_PROCESSING.md`, `code_reference/`)
- Tests: `test/test_regression.py`, `test/test_parsing.py`; inputs in
  `test/sample_off_files/`, baselines in `test/baseline_outputs/`, FF comparison fixtures
  in `test/compare/`
- Main entry: `src/main.py` (for testing/development)
- Task tracking: `TASKS.md` (current work), `archive/` (completed work periods)

### TASKS.md Archive System

**Purpose**: Archive completed TASKS.md files to maintain a clean, focused current task list while preserving historical work records in git.

**Archive location**: `archive/` directory (tracked in git)

**Filename format**: `TASKS_YYYY-MM-DD_to_YYYY-MM-DD.md`
- Example: `TASKS_2025-12-03_to_2025-12-05.md`

**Required summary line**: Each archived file must include a greppable summary at the top:
```markdown
<!-- ARCHIVE_SUMMARY: Brief description of work (YYYY-MM-DD to YYYY-MM-DD) -->
```

**Quick search** for all archived work:
```bash
grep "ARCHIVE_SUMMARY:" archive/*.md
```

**When to archive**:
- Current TASKS.md becomes too long or cluttered
- Major milestone or work period completed
- Starting a new phase of development

**How to archive**:
1. Copy TASKS.md to `archive/TASKS_START-DATE_to_END-DATE.md`
2. Add `ARCHIVE_SUMMARY` line at the top with brief description
3. Create fresh TASKS.md with empty sections
4. Commit and push to git

See `archive/README.md` for detailed archiving instructions and examples.

### Important Atom Types

The package filters out special atom types:
- "NETF": Net force atoms (excluded from most operations)
- "TORQ": Torque atoms (excluded from most operations)

Use `functions._remove_netf_torq_atname()` and `functions._remove_netf_torq_atnum()` to filter these.

### Error Handling

Functions use `exit(1)` for critical errors. When modifying code, maintain this pattern for consistency and provide clear error messages describing the issue and suggesting fixes.

### Name Translation

Many functions accept `name_translation` dictionary to map atom names from .off file to topology file names:
```python
name_translation = {'OffAtom1': 'TopAtom1', 'OffAtom2': 'TopAtom2'}
```

### Filtering Options

Both `off.gmx.gen_nonbonded_tabpot()` and `off.gmx.gen_nonbonded_topology()` support:
- `incl_mol`: List of molecule names to include (default: all)
- `excl_interactions`: List of interaction names to exclude (must match .off file exactly)
- `excl_pairs`: List of atom pairs `[['At1', 'At2'], ...]` to exclude

`off.openmm.gen_xml()` supports `incl_mol` only.

## Testing

Run the suite from the repository root:
```bash
python -m pytest test/ -v
```
Both test modules insert `src/` into `sys.path` themselves, so no install is needed.

- `test/test_regression.py`: tests 1-9. Tests 1-8 regenerate GROMACS outputs into a
  `tmp_path` and diff them against `test/baseline_outputs/testN_*/` (`.xvg` compared
  numerically with `numpy.allclose`, `.top` line-by-line). `test9_change_molecule_blypsp4f`
  is assertion-based and has no baseline directory.
- `test/test_parsing.py`: compares `off.nonbonded` / `off.bonded` against JSON baselines
  in `test/baseline_outputs/parsing/`.

See `TESTING.md` for per-test detail, the `template.top` format requirements, and how to
regenerate baselines.

Sample .off files available in `test/sample_off_files/`:
- `methane_intra.off`: Simple single-molecule test case
- `ethane_intra.off`: Small hydrocarbon test case
- `water_intra.off`: Water model test case
- `butanediol_intra.off`: Medium complexity molecule
- `big_alanine.off`: Large biomolecule test case
- `curcubituril.off`: Macrocycle test case
- `h_butanol_fitwater.off`: Solute with co-fitted water (used by `change_molecule` tests)

Force field comparison files in `test/compare/`: `base.off` and `compare.off`

**Not covered by tests**: `populate_bonded()` has no test and no `atoms.dat`/`bonds.dat`
fixtures in the repo. Anything you change there is unguarded — add a fixture if you touch it.
