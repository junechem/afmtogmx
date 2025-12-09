# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`afmtogmx` is a Python package that converts .off files produced by CRYOFF into input files for GROMACS molecular dynamics simulations. The package generates tabulated potentials and topology files for both bonded and nonbonded interactions.

## Core Architecture

### Main Entry Point

The primary class is `ReadOFF` in `src/afmtogmx/core/gen_md.py`. All workflow operations start by instantiating this class:

```python
import afmtogmx as afm
off = afm.ReadOFF(off_loc="path/to/intra.off")
```

### Key Data Structures

After initialization, `ReadOFF` objects contain:
- `off.bonded`: Dictionary with fitted bonded parameters organized by molecule name and interaction type (ATO, BON, ANG, BD3, DIH, CDI)
- `off.nonbonded`: Dictionary with nonbonded parameters organized by atom pairs and interaction types
- `off.charges`: Dictionary with atomic charges per molecule (default 0.0 for all atoms; must be manually set if needed)
- `off.residues`: Dictionary with residue definitions and atom number mappings
- `off.sections`: Internal dictionary splitting .off file into 5 sections (ff_input, intra_potential, inter_potential, molecular_definition, table_potential)
- `off.config`: Dictionary storing default parameter values for all workflow methods
- `off.nonbonded_tabpot`: Automatically populated by `gen_nonbonded_tabpot()` (None until generated)
- `off.bonded_tabpot`: Automatically populated by `gen_bonded_tabpot()` (None until generated)

### Configuration System

The `ReadOFF` class includes a configuration system that allows setting default parameters once instead of passing them to each method call.

**Setting configuration values**:
```python
off.set_config(spacing_nonbonded=0.001, scale_C6=False, tabpot_prefix='my_table')
# Returns self for method chaining
```

**Getting configuration values**:
```python
all_config = off.get_config()  # Returns full config dictionary
spacing = off.get_config('spacing_nonbonded')  # Returns specific value
```

**Configuration keys**:
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
- `name_translation`: Atom name translation dictionary (default: {})
- `special_pairs`: Custom attractive interactions per pair (default: empty dict)

**Parameter resolution order**:
1. Explicitly passed parameter (if not None)
2. Configuration value
3. Built-in default

This means you can set defaults via `set_config()` and override them on individual method calls as needed.

### Module Responsibilities

**gen_md.py** (495 lines)
- Main `ReadOFF` class with all public API methods
- Parses .off files and populates data structures
- Orchestrates calls to other modules

**functions.py** (550 lines)
- Low-level parsing utilities for .off file format
- Interaction function definitions (exp, srd, shtr, powe, quarbond)
- Helper functions for filtering and organizing parameters
- Maintains global counter `total_bonded_added` for tracking

**tabulated_potentials.py** (323 lines)
- Generates tabulated potential tables for bonded and nonbonded interactions
- Handles special interaction types: POW_*, DPO_*, SRD_*, PEX_*, EXP, BUC, etc.
- C6 scaling for dispersion corrections
- Soft-core transformations for free energy calculations

**topology.py** (448 lines)
- Generates GROMACS .top topology files
- Writes [ nonbonded_params ] and [ pairs ] sections
- Template-based topology file modification
- Atom name translation support

**residues.py** (121 lines)
- Manages residue definitions and atom number mappings
- Validates user-provided residue definitions

**compare.py** (67 lines)
- Compares two force fields (two ReadOFF objects)
- Generates difference strings for documentation

## Standard Workflow

The typical workflow follows this sequence:

**Traditional approach** (passing parameters to each method):

1. **Read .off file**: `off = afm.ReadOFF(off_loc="path/to/intra.off")`
2. **Set charges manually** (if needed): `off.charges['MolName']['AtomName'] = charge_value`
3. **Generate tabulated potentials**:
   - `nonbonded_tabpot = off.gen_nonbonded_tabpot(...)`
   - `bonded_tabpot = off.gen_bonded_tabpot(...)` (if needed)
4. **Write tabulated potentials**:
   - `off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot)`
   - `off.write_bonded_tabpot(bonded_tabpot=bonded_tabpot)`
5. **Generate topology files**:
   - `off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')`
   - `off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')`

**Configuration-based approach** (recommended for cleaner code):

1. **Read .off file and configure**:
   ```python
   off = afm.ReadOFF(off_loc="path/to/intra.off")
   off.set_config(tabpot_prefix='my_table', tabpot_dir='tables/', scale_C6=False)
   ```
2. **Load charges**: `off.load_charges_from_file('charges.txt')`
3. **Generate and write tabulated potentials**:
   ```python
   # Generate methods store results automatically in off.nonbonded_tabpot and off.bonded_tabpot
   off.gen_nonbonded_tabpot()  # Uses config values, stores in off.nonbonded_tabpot
   off.gen_bonded_tabpot()     # Uses config values, stores in off.bonded_tabpot

   # Write methods use stored results automatically - no variables needed!
   off.write_nonbonded_tabpot()  # Uses off.nonbonded_tabpot and config prefix/dir
   off.write_bonded_tabpot()     # Uses off.bonded_tabpot and config prefix/dir
   ```
4. **Generate topology files**:
   ```python
   off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')
   off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')
   ```

**Note**: You can still pass custom dictionaries to write methods for advanced use cases:
```python
custom_tabpot = off.gen_nonbonded_tabpot()  # Can still capture if needed
off.write_nonbonded_tabpot(nonbonded_tabpot=custom_tabpot)  # Can still pass explicitly
```

See detailed workflow examples in `src/afmtogmx/__init__.py`.

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

### Soft-Core Scaling

Use `sc_sigma` parameter in `gen_nonbonded_tabpot()` for free energy calculations. Tables are scaled to work correctly with the sc-sigma value in GROMACS .mdp files.

### Recent Changes

**Configuration System (Latest)**:
- Added configuration system with `set_config()` and `get_config()` methods
- All workflow methods now support parameter resolution (explicit → config → default)
- Added `load_charges_from_file()` convenience method for loading charges from text files
- Converted `write_nonbonded_tabpot()` and `write_bonded_tabpot()` from static methods to instance methods

**Breaking Changes**:
- `gen_nonbonded_tabpot()`: Parameter `spacing` renamed to `spacing_nonbonded`, `length` renamed to `length_nonbonded`
- `gen_bonded_tabpot()`: Parameter `spacing` renamed to `spacing_bonded`, `length` renamed to `length_bonded`
- `write_nonbonded_tabpot()` and `write_bonded_tabpot()`: Parameters `prefix` renamed to `tabpot_prefix`, `to_dir` renamed to `tabpot_dir`
- All workflow methods now use `None` as default for most parameters to enable config resolution

**Previous Changes**:
- Bug fix for `gen_nonbond_tabpam` function (commit 6289640)
- Added C12 column scaling for HFE (Hydration Free Energy) calculations (commits 5c102af, 96f8e83)
- Added `excl_pairs` functionality to exclude specific atom pairs from tabulated potentials
- Removed all charge calculation functions - only manual charge assignment is now supported

## Development Notes

### File Organization

- Core logic: `src/afmtogmx/core/`
- Package exports: `src/afmtogmx/__init__.py`
- Test files: `test/sample_off_files/` (sample .off files), `test/compare/` (force field comparison)
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

Both `gen_nonbonded_tabpot()` and `gen_nonbonded_topology()` support:
- `incl_mol`: List of molecule names to include (default: all)
- `excl_interactions`: List of interaction names to exclude (must match .off file exactly)
- `excl_pairs`: List of atom pairs `[['At1', 'At2'], ...]` to exclude

## Testing

Sample .off files available in `test/sample_off_files/`:
- `methane_intra.off`: Simple single-molecule test case
- `ethane_intra.off`: Small hydrocarbon test case
- `water_intra.off`: Water model test case
- `butanediol_intra.off`: Medium complexity molecule
- `big_alanine.off`: Large biomolecule test case

Force field comparison files in `test/compare/`: `base.off` and `compare.off`
