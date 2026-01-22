# afmtogmx

A Python package that converts force field parameters from [CRYOFF](https://wanglab.hosted.uark.edu/cryoff/wanglab_CRYOFF.html) `.off` files into input files for [GROMACS](https://www.gromacs.org/) molecular dynamics simulations.

## Overview

**afmtogmx** bridges the gap between quantum-mechanically derived force fields and production molecular dynamics simulations. It takes fitted interaction parameters from CRYOFF—a force field optimization tool that fits parameters to ab initio data—and generates the tabulated potentials and topology files required by GROMACS.

### What It Does

- Parses CRYOFF `.off` files containing bonded and nonbonded interaction parameters
- Generates tabulated potential files (`.xvg`) for custom interaction functions (exponential, power-law, etc.)
- Creates GROMACS topology files (`.top`)
- Supports features like C6 dispersion scaling, soft-core transformations for free energy calculations, and custom atom name translations

## Installation

Clone the repository and install in development mode:

```bash
git clone https://github.com/junechem/afmtogmx.git
cd afmtogmx
pip install -e .
```

### Requirements

- Python 3.8+
- NumPy 1.20+

## Quick Start

```python
import afmtogmx

# 1. Load the .off file
off = afmtogmx.ReadOFF(off_loc="path/to/forcefield.off")

# 2. Configure settings (optional)
off.set_config(
    tabpot_dir='tables/',
    tabpot_prefix='table',
    scale_C6=True
)

# 3. Load charges from file (if needed)
off.load_charges_from_file('charges.txt')

# 4. Generate and write tabulated potentials
off.gen_nonbonded_tabpot()
off.gen_bonded_tabpot()
off.write_nonbonded_tabpot()
off.write_bonded_tabpot()

# 5. Generate topology files
off.gen_nonbonded_topology(template_file='template.top', write_to='temp.top')
off.gen_bonded_topology(template_file='temp.top', write_to='topol.top')
```

## Usage Guide

### Reading Force Field Data

After loading an `.off` file, parameters are accessible as dictionaries:

```python
off = afmtogmx.ReadOFF(off_loc="methane.off")

# Access bonded parameters (organized by molecule → interaction type)
print(off.bonded.keys())           # e.g., dict_keys(['CH4'])
print(off.bonded['CH4'].keys())    # e.g., dict_keys(['ATO', 'BON', 'ANG', ...])

# Access nonbonded parameters (organized by atom pair → interaction type)
print(off.nonbonded.keys())        # e.g., dict_keys([('C', 'H'), ('H', 'H'), ...])
```

### Setting Charges

Charges default to zero and must be set manually:

```python
# Option 1: Direct assignment
off.charges['H2O']['O'] = -0.82
off.charges['H2O']['H'] = 0.41

# Option 2: Load from file
off.load_charges_from_file('charges.txt')
```

Charge file format:
```
# Comment lines start with #
H2O
O -0.82
H 0.41
```

### Configuration System

Set default parameters once instead of passing them to each method:

```python
off.set_config(
    spacing_nonbonded=0.0005,   # Table spacing in nm
    length_nonbonded=3.0,       # Table length in nm
    scale_C6=True,              # Scale for dispersion corrections
    tabpot_prefix='table',      # Output file prefix
    tabpot_dir='potentials/'    # Output directory
)

# Override individual settings when needed
off.gen_nonbonded_tabpot(spacing_nonbonded=0.001)
```

### Filtering Interactions

Control which interactions are included in the output:

```python
off.gen_nonbonded_tabpot(
    incl_mol=['MOL1', 'MOL2'],           # Only include specific molecules
    excl_interactions=['BUCWATER'],       # Exclude specific interaction types
    excl_pairs=[['O', 'O'], ['H', 'H']]  # Exclude specific atom pairs
)
```

### Atom Name Translation

Map atom names from the `.off` file to topology names:

```python
off.set_config(name_translation={
    'C_off': 'C_gro',
    'H_off': 'H_gro'
})
```

## API Reference

### `ReadOFF(off_loc)`

Main class for parsing `.off` files and generating GROMACS inputs.

**Attributes:**
| Attribute | Description |
|-----------|-------------|
| `bonded` | Nested dict of bonded parameters by molecule and interaction type |
| `nonbonded` | Dict of nonbonded parameters by atom pair |
| `charges` | Dict of atomic charges by molecule and atom name |
| `config` | Configuration settings dictionary |
| `nonbonded_tabpot` | Generated nonbonded potentials (after calling `gen_nonbonded_tabpot()`) |
| `bonded_tabpot` | Generated bonded potentials (after calling `gen_bonded_tabpot()`) |

**Methods:**

| Method | Description |
|--------|-------------|
| `set_config(**kwargs)` | Set configuration values; returns `self` for chaining |
| `get_config(key=None)` | Get all config or a specific value |
| `load_charges_from_file(filepath)` | Load charges from a text file |
| `gen_nonbonded_tabpot(...)` | Generate nonbonded tabulated potentials |
| `gen_bonded_tabpot(...)` | Generate bonded tabulated potentials |
| `write_nonbonded_tabpot(...)` | Write nonbonded `.xvg` files |
| `write_bonded_tabpot(...)` | Write bonded `.xvg` files |
| `gen_nonbonded_topology(...)` | Generate topology with nonbonded parameters |
| `gen_bonded_topology(...)` | Generate topology with bonded parameters |

### Configuration Options

| Key | Default | Description |
|-----|---------|-------------|
| `spacing_nonbonded` | 0.0005 | Table spacing for nonbonded potentials (nm) |
| `length_nonbonded` | 3.0 | Table length for nonbonded potentials (nm) |
| `spacing_bonded` | 0.0001 | Table spacing for bonded potentials (nm) |
| `length_bonded` | 0.3 | Table length for bonded potentials (nm) |
| `scale_C6` | True | Scale columns 4-5 by C6 for dispersion corrections |
| `sc_sigma` | 0.0 | Soft-core sigma for free energy calculations |
| `tabpot_prefix` | 'table' | Prefix for output potential files |
| `tabpot_dir` | '' | Directory for output files |
| `incl_mol` | [] | Molecules to include (empty = all) |
| `excl_interactions` | [] | Interaction types to exclude |
| `excl_pairs` | [] | Atom pairs to exclude |
| `name_translation` | {} | Atom name mapping dictionary |
| `special_pairs` | {} | Custom attractive interactions per pair |

## Supported Interaction Types

**Nonbonded:**
- `EXP` — Exponential repulsion
- `BUC` — Buckingham potential
- `POW_n` — Power-law (e.g., `POW_6`, `POW_12`)
- `DPO_n` — Damped power-law
- `SRD_n` — Short-range damped
- `PEX_n` — Power-exponential
- `COU` — Coulombic

**Bonded:**
- `BON` — Bond stretching
- `ANG` — Angle bending
- `DIH` — Dihedral torsion
- `BD3` — Out-of-plane bending
- `CDI` — Coupled dihedral

## Project Structure

```
afmtogmx/
├── src/afmtogmx/
│   ├── __init__.py          # Package exports and documentation
│   └── core/
│       ├── gen_md.py        # Main ReadOFF class
│       ├── functions.py     # Parsing utilities
│       ├── tabulated_potentials.py  # Potential generation
│       ├── topology.py      # Topology file generation
│       └── residues.py      # Residue management
├── test/
│   ├── sample_off_files/    # Example .off files
│   └── baseline_outputs/    # Test cases
└── pyproject.toml
```

## License

MIT License

## Related Projects

- [CRYOFF](https://wanglab.hosted.uark.edu/cryoff/wanglab_CRYOFF.html) — Force field optimization from ab initio data
- [GROMACS](https://www.gromacs.org/) — Molecular dynamics simulation package