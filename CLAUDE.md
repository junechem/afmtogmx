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
- `off.charges`: Dictionary with atomic charges per molecule (default 0.0 until `calc_charges()` is called)
- `off.residues`: Dictionary with residue definitions and atom number mappings
- `off.sections`: Internal dictionary splitting .off file into 5 sections (ff_input, intra_potential, inter_potential, molecular_definition, table_potential)

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

**chargefxns.py** (192 lines)
- Derives atomic charges from coulombic interactions in .off file
- Charge normalization methods (currently "M-POPULOUS")
- Residue-based charge neutralization
- Works from known atom charges or derives from self-interactions

**residues.py** (121 lines)
- Manages residue definitions and atom number mappings
- Validates user-provided residue definitions
- Generates residue priority for charge neutralization

**compare.py** (67 lines)
- Compares two force fields (two ReadOFF objects)
- Generates difference strings for documentation

## Standard Workflow

The typical workflow follows this sequence:

1. **Read .off file**: `off = afm.ReadOFF(off_loc="path/to/intra.off")`
2. **Calculate charges** (if coulombic interactions present): `off.calc_charges(...)`
3. **Generate tabulated potentials**:
   - `nonbonded_tabpot = off.gen_nonbonded_tabpot(...)`
   - `bonded_tabpot = off.gen_bonded_tabpot(...)` (if needed)
4. **Write tabulated potentials**:
   - `off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot)`
   - `off.write_bonded_tabpot(bonded_tabpot=bonded_tabpot)`
5. **Generate topology files**:
   - `off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')`
   - `off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')`

See detailed workflow examples in `src/afmtogmx/__init__.py` (lines 58-68).

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

### Charge Calculation

Charges are derived from coulombic (COU) interactions in the .off file using:
- A known atom with known charge, OR
- A known atom with self-interaction and charge sign (derives charge as ±√(self-interaction))

Charges are normalized using "M-POPULOUS" method: excess charge distributed to most populous atom type with nonzero charge.

### Residue-Based Neutralization

Set `neutral_residues=True` in `calc_charges()` to neutralize individual residues rather than entire molecules. Use `residue_priority` dictionary to control which residues are neutralized first when atom types are shared.

### Soft-Core Scaling

Use `sc_sigma` parameter in `gen_nonbonded_tabpot()` for free energy calculations. Tables are scaled to work correctly with the sc-sigma value in GROMACS .mdp files.

### Recent Changes

Based on recent commits:
- Bug fix for `gen_nonbond_tabpam` function (commit 6289640)
- Added C12 column scaling for HFE (Hydration Free Energy) calculations (commits 5c102af, 96f8e83)
- Added `excl_pairs` functionality to exclude specific atom pairs from tabulated potentials
- Charge normalization rewrite in progress (`chargefxns_old.py` is previous version)

## Development Notes

### File Organization

- Core logic: `src/afmtogmx/core/`
- Package exports: `src/afmtogmx/__init__.py`
- Test files: `test/sample_off_files/` (sample .off files), `test/compare/` (force field comparison)
- Main entry: `src/main.py` (for testing/development)

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
