# afmtogmx

A Python package that converts force field parameters from [CRYOFF](https://wanglab.hosted.uark.edu/cryoff/wanglab_CRYOFF.html) `.off` files into simulation-ready input files for [GROMACS](https://www.gromacs.org/) and [OpenMM](https://openmm.org/).

## Overview

**afmtogmx** bridges the gap between quantum-mechanically derived force fields and production molecular dynamics simulations. It takes fitted interaction parameters from CRYOFF—a force field optimization tool that fits parameters to ab initio data—and generates the input files each engine needs.

### What It Does

- Parses CRYOFF `.off` files containing bonded and nonbonded interaction parameters
- **GROMACS**: generates tabulated potential files (`.xvg`) for custom interaction functions and writes complete topology files (`.top`)
- **OpenMM**: generates a self-contained `<ForceField>` XML and preprocesses PDB files to match it
- Rewrites the force field before output: swap a molecule for a stored reference force field, or rebuild a molecule's bonded terms from an atom/bond topology
- Exports a publication-ready text summary of every fitted parameter
- Supports C6 dispersion scaling, soft-core transformations for free energy calculations, and atom name translation

## Installation

Clone the repository and install in development mode:

```bash
git clone https://github.com/junechem/afmtogmx.git
cd afmtogmx
pip install -e ".[dev]"     # omit [dev] if you don't need pytest
```

A conda environment file is also provided:

```bash
conda env create -f environment.yml
conda activate afmtogmx
```

### Requirements

- Python 3.8+
- NumPy 1.20+
- OpenMM — **optional**. It is only used to look up atomic masses when writing the
  force field XML. Without it, XML generation still works but every `mass` attribute
  is written as `0.0`, which you must fill in before running a simulation.

## Architecture at a Glance

Parsing a `.off` file gives you a `ReadOFF` object. It holds the force field data, and
exposes one backend per simulation engine:

```
off = afmtogmx.ReadOFF('forcefield.off')

  off.bonded / off.nonbonded / off.charges / off.residues   ← shared parsed data
    │
    ├── off.gmx      → tabulated potentials (.xvg) + topology (.top)
    └── off.openmm   → force field XML + processed PDB
```

Anything that modifies the force field itself (charges, `change_molecule`,
`populate_bonded`, `gen_residues`) lives directly on `off` and affects both backends.
Anything engine-specific lives on the backend.

## Quick Start — GROMACS

```python
import afmtogmx

off = afmtogmx.ReadOFF(off_loc='ethane_intra.off')
off.load_charges_from_file('charges.txt')

off.gmx.set_config(tabpot_dir='tables/', tabpot_prefix='table', scale_C6=True)

# Generated tables are stored on the backend, so nothing needs passing around
off.gmx.gen_nonbonded_tabpot()      # → off.gmx.nonbonded_tabpot
off.gmx.gen_bonded_tabpot()         # → off.gmx.bonded_tabpot
off.gmx.write_nonbonded_tabpot()
off.gmx.write_bonded_tabpot()

# Topology is built in two passes over a template
off.gmx.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')
off.gmx.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')
```

`write_nonbonded_tabpot()` prints the unique atom list and the pair list, which are what
you paste into your `.mdp` as `energygrps` and `energygrp_table`.

The `template.top` you pass in must already contain the `[ atoms ]`, `[ bonds ]`,
`[ angles ]`, `[ dihedrals ]` and `[ exclusions ]` section headers, separated by blank
lines — afmtogmx fills them in, it does not create them. See [TESTING.md](TESTING.md) for
a worked template.

## Quick Start — OpenMM

```python
import afmtogmx

off = afmtogmx.ReadOFF(off_loc='butanol_water.off')
off.load_charges_from_file('charges.txt')

# Map .off molecule names onto the residue names used in your PDB
off.openmm.set_config(molname_translations={'H2OQM': 'SOL'})
off.openmm.gen_xml(output='forcefield.xml')

# Rename PDB atoms to match the XML and emit fresh CONECT records
off.openmm.preprocess_pdb('conf.pdb', 'forcefield.xml', 'conf_processed.pdb')
```

`forcefield.xml` and `conf_processed.pdb` load directly into an OpenMM simulation — the
XML is self-contained, so there is no second force field file to merge and no naming
conflicts to resolve.

A few things worth knowing:

- Atom types are namespaced `<MOLNAME>_<TYPE>` (e.g. `H2OQM_OW`) to avoid collisions
  between molecules.
- Units are converted for you: kcal/mol → kJ/mol, Å → nm.
- The XML path supports a narrower set of nonbonded interactions than the GROMACS path:
  `EXP`, `STR`/`STRC`, `SRD`, `POW` (emitted as `SRD` with `r0 = 0`), and `BUC` (split
  into `EXP` + `SRD`).
- `preprocess_pdb` reads **only** the XML, never in-memory state, so one generated XML
  can process any number of PDB files. It uses connectivity (graph isomorphism) matching
  when the input PDB has CONECT records — which also reorders atoms into the XML's
  canonical order — and falls back to positional matching otherwise, guarded by a
  `maxwarn` element-mismatch failsafe.

Full detail on the PDB stage is in [docs/PDB_PROCESSING.md](docs/PDB_PROCESSING.md).

## Setting Charges

Charges default to zero and must be set explicitly:

```python
# Option 1: direct assignment
off.charges['H2O']['O'] = -0.82
off.charges['H2O']['H'] = 0.41

# Option 2: load from a file (chainable)
off.load_charges_from_file('charges.txt')
```

The charge file lists a molecule name on its own line, then one `atom charge` pair per
line. Blank lines and `#` comments are ignored:

```
# TIP3P-like water
H2O
O -0.82
H 0.41
```

A legacy format with no molecule names is also accepted — a bare `atom charge` line
appearing before any molecule name is applied to *every* molecule containing that atom:

```
OW -0.82
HW  0.41
```

Atoms not named in the file keep their default charge of 0.0. Atoms that are named are
overwritten.

## Configuring a Backend

Each backend keeps its own config, so you set defaults once instead of repeating them on
every call:

```python
off.gmx.set_config(spacing_nonbonded=0.001, scale_C6=False, tabpot_dir='tables/')
off.gmx.get_config()              # full dict
off.gmx.get_config('scale_C6')    # one value

# Config values can still be overridden per call
off.gmx.gen_nonbonded_tabpot(spacing_nonbonded=0.0005)
```

Resolution order is: explicitly passed argument → config value → built-in default.
`set_config()` returns the backend, so calls chain.

### `off.gmx` config keys

| Key | Default | Description |
|-----|---------|-------------|
| `spacing_nonbonded` | `0.0005` | Table spacing for nonbonded potentials (nm) |
| `length_nonbonded` | `3` | Table length for nonbonded potentials (nm) |
| `spacing_bonded` | `0.0001` | Table spacing for bonded potentials (nm) |
| `length_bonded` | `0.3` | Table length for bonded potentials (nm) |
| `scale_C6` | `True` | Scale columns 4–5 by C6 so GROMACS dispersion corrections work |
| `sc_sigma` | `0.0` | Soft-core sigma for free energy calculations |
| `tabpot_prefix` | `'table'` | Prefix for output `.xvg` files |
| `tabpot_dir` | `''` | Directory for output `.xvg` files |
| `write_blank` | `True` | Also write the blank (all-zero) table |
| `incl_mol` | `[]` | Molecules to include (empty = all) |
| `excl_interactions` | `[]` | Interaction types to exclude |
| `excl_pairs` | `[]` | Atom pairs to exclude |
| `name_translation` | `{}` | Atom name mapping, `.off` name → topology name |
| `special_pairs` | `{}` | Per-pair override of which interactions count as attractive |

By default the attractive interactions (`POW_6`, `DPO_6`, `SRD_6`) go into columns 4–5 of
each table and everything except `COU` goes into columns 6–7. `special_pairs` overrides
that per pair — e.g. `{('At1', 'At2'): ['POW_6', 'POW_8']}` — and is incompatible with
`scale_C6=True`, which assumes exactly one attractive interaction per pair.

### `off.openmm` config keys

| Key | Default | Description |
|-----|---------|-------------|
| `incl_mol` | `[]` | Molecules to include (empty = all) |
| `molname_translations` | `{}` | Molecule name mapping, `.off` name → PDB residue name |

> **Three similar names, three different things.** `name_translation` (`off.gmx` config)
> maps *atom* names. `molname_translations` (`off.openmm` config, plural) maps *molecule*
> names to PDB residue names. `molname_translation` (singular) is a `write_report()`
> keyword argument for display names in the report. They are not interchangeable.

## Filtering What Gets Written

Both backends accept the same filters, either in config or per call:

```python
off.gmx.gen_nonbonded_tabpot(
    incl_mol=['MOL1', 'MOL2'],            # only these molecules
    excl_interactions=['BUCWATER'],        # names must match the .off file exactly
    excl_pairs=[['O', 'O'], ['H', 'H']],   # skip these atom pairs
)

off.openmm.gen_xml(incl_mol=['H2OQM'], output='water_only.xml')
```

Atom name translation maps `.off` names onto the names your topology uses:

```python
off.gmx.set_config(name_translation={'C_off': 'C_gro', 'H_off': 'H_gro'})
```

## Replacing a Molecule with a Reference Force Field

An AFM-fitted `.off` file often contains refitted solvent parameters you don't want to
keep. `change_molecule()` swaps that molecule's parameters for a stored reference force
field, in place, before any output is generated:

```python
off = afmtogmx.ReadOFF(off_loc='h_butanol_fitwater.off')

off.change_molecule(
    mol_name='H2OQM',            # molecule as named in the loaded .off file
    reference_ff='BLYPSP-4F',    # file stem in src/afmtogmx/reference_ff/
    ref_mol_name='H2OQM',        # which molecule to take from the reference
)
```

What it does:

- Replaces the molecule's bonded parameters wholesale.
- Replaces the pure solvent–solvent nonbonded pairs with the reference's.
- **Leaves cross-term pairs alone.** Solute–solvent parameters were fitted against the
  reference water model, so they keep their original fitted values; only the atom type
  name on the replaced molecule's side is renamed.
- Auto-loads a sibling `.charges` file if the reference ships one.
- Rebuilds `off.residues`, and returns `self` so the call chains.

`ref_mol_name` is only optional when the reference force field defines exactly one
molecule. The bundled `BLYPSP-4F` defines several, so it is required there.

Pass `atom_name_map` to rename types at the same time — useful when the reference's names
would otherwise collide with your solute's:

```python
off.change_molecule(
    mol_name='H2OQM',
    reference_ff='BLYPSP-4F',
    ref_mol_name='H2OQM',
    atom_name_map={'OW': 'OW_sp4f', 'HW': 'HW_sp4f'},
)
```

Reference force fields live in `src/afmtogmx/reference_ff/`. `BLYPSP-4F` ships with the
package; add more by dropping a `<name>.off` (and optional `<name>.charges`) there.

## Building a Molecule from Atoms and Bonds

CRYOFF often fits a small monomer, but you want to simulate a large assembled molecule
built from the same atom types. `populate_bonded()` rebuilds a molecule's entire bonded
section from a topology you supply as atoms and bonds — angles, dihedrals and exclusions
are derived from the bond graph:

```python
off = afmtogmx.ReadOFF(off_loc='cucurbituril/intra.off')   # monomer + co-fitted water

new_off = off.populate_bonded(
    directory='cb7_topology/',    # holds the four files described below
    new_mol_name='CB7',           # name for the assembled molecule
    remove_molecules=['C7F'],     # drop the original monomer's entries
)

new_off.load_charges_from_file('cb7_charges.txt')   # new molecule starts at 0.0
new_off.gmx.gen_nonbonded_tabpot()
```

This returns a **new** `ReadOFF`; the original `off` is left untouched. Co-fitted
molecules (water, here) carry through unchanged, and `nonbonded` is pruned to pairs whose
atom types still appear in a surviving molecule.

### Input directory

Four files, all required:

| File | Format |
|------|--------|
| `atoms.dat` | `atom_number  atom_type` |
| `bonds.dat` | `atom_a  atom_b` |
| `valid_dihedrals.dat` | `type1-type2-type3-type4`, one signature per line |
| `parameters.dat` | `#define <signature> <subtype> <values...>` |

`parameters.dat` looks like this — the number of underscores in the signature selects the
term kind (1 = bond, 2 = angle, 3 = dihedral):

```
#define   C1_C1         HAR   0.15468521   149286.94
#define   C1_C1_H1      HAR   110.04697    384.75972
#define   C1_C1_N1_C2   NCO   0.0          -8.137997   3
```

Things that will otherwise bite you:

- **`parameters.dat` values are in GROMACS units** (nm, kJ/mol, dihedral order
  `phi K mult`). Conversion to the `.off`'s internal Å/kcal units is automatic.
- **`valid_dihedrals.dat` is a filter.** Dihedrals found in the bond graph whose signature
  isn't listed are silently dropped — that is the file's purpose. Conversely, every
  signature you *do* list must resolve to a `#define`, even if no matching dihedral exists.
- Supported subtypes: bonds `HAR`; angles `HAR`; dihedrals `NCO`, `COS`, `HAR`.
- Exclusions are generated for 1-2 and 1-3 neighbours.
- Improper/out-of-plane (`BD3`) and coupled dihedral (`CDI`) terms are **not** generated.
- The call **fails fast**: every missing file, unknown atom type, and unresolved
  bond/angle/dihedral signature is collected into a single `ValueError`, so you fix all
  the gaps in one pass rather than one error at a time.

## Defining Residues

By default each molecule is one residue. `gen_residues()` splits a molecule into named
residues, either by atom type or by explicit atom numbers:

```python
off.gen_residues(
    residue_definition={'PROT': {'ALA': ['CA', 'CB', 'HA'], 'GLY': ['CA', 'HA']}},
)

off.gen_residues(
    residue_atnums={'PROT': {'ALA': [[1, 2, 3, 4], [5, 6, 7, 8]]}},
)
```

Definitions are validated against the parsed `off.bonded` data, so a typo in an atom type
is caught immediately.

## Exporting a Parameter Report

`write_report()` writes a fixed-width text summary of charges, bonded parameters, and
nonbonded parameters — sized to paste into a Word document at default margins with
Courier New 10 pt:

```python
off.write_report(
    'forcefield_parameters.txt',
    incl_mol=['UNK'],
    name_translation={'C1': 'C_alpha'},        # display names for atoms
    molname_translation={'UNK': 'butanol'},    # display names for molecules
    notation='standard',
)
```

Units are the `.off` file's native kcal/Å — nothing is converted. `notation='standard'`
uses publication-style column labels (`A`, `alpha`, `Cn`, `R0`); `notation='PN'` uses the
CRYOFF manual's generic positional labels (`P1`, `P2`, …). Nonbonded pairs appear under a
molecule's section if *either* atom belongs to it, so solute–solvent cross terms are not
lost when you filter with `incl_mol`.

## Reading Force Field Data

Everything the parser produces is a plain dictionary, so you can inspect or edit it
directly:

```python
off = afmtogmx.ReadOFF(off_loc='methane_intra.off')

off.bonded.keys()                # molecule names, e.g. dict_keys(['UNK'])
off.bonded['UNK'].keys()         # dict_keys(['ATO', 'BON', 'ANG', 'BD3', 'DIH', 'CDI', 'EXC'])
off.nonbonded.keys()             # atom pairs, e.g. dict_keys([('H', 'H'), ('C', 'C'), ('C', 'H')])
off.charges['UNK']               # {'C': 0.0, 'H': 0.0}
off.residues                     # {'Definitions': {...}, 'Residues': {...}}
```

## Supported Interaction Types

**Nonbonded (van der Waals)**: `GLJ`, `BUC`, `DBU`, `STR`, `EXP`, `GEX`, `POW`, `CSP`,
`GDP`, `PEX`, `DPO`, `SRD`, `TTP`

**Nonbonded (electrostatic)**: `COU`, `THC`

**Bonded**: bonds `HAR`, `QUA`; angles `HAR`, `QUA`, `RAH`; three-body cross terms
`QBB`, `MUB`; dihedrals `HAR`, `NCO`, `COS`; coupled dihedrals `CNCO`, `CCOS`

Two conventions to be aware of when reading `off.nonbonded` keys:

- Only the **first three characters** name the functional form. Anything after that is a
  user-chosen identifier from whoever wrote the fit — `EXPW`, `EXPINTER` and `EXPINTRA`
  are all `EXP`, and `STRC` is `STR`.
- Variable-power types are stored internally as `TYPE_POWER`, e.g. `POW_6`, `SRD_8`,
  `DPO_6`.

[docs/CRYOFF_REFERENCE.md](docs/CRYOFF_REFERENCE.md) is the authoritative source for what
every symbol means, its functional form, parameter order, and units — including the
naming inconsistencies this package carries (`CDI` vs `CDIH`, `STR` vs `STRC`, `POW` vs
`POW_6`). Consult it before adding or changing an interaction type.

## Project Structure

```
afmtogmx/
├── src/afmtogmx/
│   ├── __init__.py
│   ├── core/
│   │   ├── gen_md.py               # ReadOFF: parsing + shared force field methods
│   │   ├── functions.py            # .off parsing utilities + interaction math
│   │   ├── gmx_backend.py          # GROMACSBackend  (off.gmx)
│   │   ├── openmm_backend.py       # OpenMMBackend   (off.openmm)
│   │   ├── tabulated_potentials.py # .xvg table generation
│   │   ├── topology.py             # .top topology generation
│   │   ├── xml_generation.py       # OpenMM <ForceField> XML sections
│   │   ├── pdb_processing.py       # PDB atom renaming + CONECT records
│   │   ├── populate_bonded.py      # bonded terms from an atoms/bonds topology
│   │   ├── report.py               # write_report text export
│   │   ├── residues.py             # residue definitions
│   │   └── compare.py              # diff two force fields
│   └── reference_ff/               # bundled reference force fields (BLYPSP-4F)
├── docs/
│   ├── CRYOFF_REFERENCE.md
│   ├── PDB_PROCESSING.md
│   └── code_reference/             # per-function reference for every core module
├── test/
│   ├── test_regression.py
│   ├── test_parsing.py
│   ├── sample_off_files/
│   ├── baseline_outputs/
│   └── compare/
└── pyproject.toml
```

## Testing

From the repository root:

```bash
python -m pytest test/ -v
```

The test modules add `src/` to `sys.path` themselves, so this works without installing
the package. `test_regression.py` regenerates the GROMACS outputs for each sample system
and diffs them against committed baselines (`.xvg` compared numerically, `.top`
line-by-line); `test_parsing.py` checks the parser against JSON baselines.

See [TESTING.md](TESTING.md) for the individual test cases, the `template.top` format
requirements, and how to regenerate baselines after an intentional change.

## Documentation

| Document | What it covers |
|----------|----------------|
| [docs/CRYOFF_REFERENCE.md](docs/CRYOFF_REFERENCE.md) | Authoritative meanings, parameter orders and units for every CRYOFF `.off` symbol |
| [docs/PDB_PROCESSING.md](docs/PDB_PROCESSING.md) | The PDB preprocessing pipeline in depth: matching strategies, edge cases, limitations |
| [docs/code_reference/](docs/code_reference/README.md) | Per-function reference for every module in `core/` |
| [TESTING.md](TESTING.md) | Test suite layout, running tests, regenerating baselines |
| [CLAUDE.md](CLAUDE.md) | Guidance for working on this codebase |

## License

MIT License

## Related Projects

- [CRYOFF](https://wanglab.hosted.uark.edu/cryoff/wanglab_CRYOFF.html) — Force field optimization from ab initio data
- [GROMACS](https://www.gromacs.org/) — Molecular dynamics simulation package
- [OpenMM](https://openmm.org/) — High-performance molecular simulation toolkit
