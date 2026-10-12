# Testing Guide

This document explains the testing infrastructure for the `afmtogmx` project.

## Running the Tests

From the repository root:

```bash
python -m pytest test/ -v            # everything
python -m pytest test/test_regression.py -v
python -m pytest test/test_parsing.py -v
python -m pytest test/test_regression.py::test9_change_molecule_blypsp4f -v
```

Both test modules insert `src/` into `sys.path` themselves, so the package does not need
to be installed. `pytest` comes from the `dev` extra: `pip install -e ".[dev]"`.

## Test Layout

```
test/
├── test_regression.py        # tests 1-9 (GROMACS outputs + change_molecule)
├── test_parsing.py           # parser output vs JSON baselines
├── sample_off_files/         # input .off files (do not modify)
├── baseline_outputs/
│   ├── parsing/              # JSON baselines for test_parsing.py
│   └── testN_*/              # per-test baselines for test_regression.py
└── compare/                  # base.off / compare.off, fixtures for core/compare.py
```

**`test_regression.py`** copies each test's inputs into a `tmp_path`, runs the workflow
there, and diffs the results against the committed baseline: `.xvg` files are compared
numerically (`numpy.allclose`, `rtol=1e-6`, `atol=1e-10`), `.top` files line-by-line with
trailing whitespace stripped. `template.top` and `charges.txt` are treated as inputs and
skipped.

**`test_parsing.py`** compares `off.nonbonded` and `off.bonded` against JSON baselines in
`test/baseline_outputs/parsing/`, parametrized over six sample `.off` files
(`methane_intra`, `ethane_intra`, `water_intra`, `butanediol_intra`, `big_alanine`,
`curcubituril`). Tuple keys are serialized as their `repr`.

## Test Structure

Each of the test 1-8 baseline directories contains:
- `generate_testX.py` - Python script to regenerate outputs
- `template.top` - Input GROMACS topology template
- `topol.top` - Final topology file (baseline)
- `temp_nonbonded.top` - Intermediate topology file (baseline)
- `tabpot/*.xvg` - Tabulated potential files (baseline)

Test 9 is assertion-based and has no baseline directory.

> **Note**: the `generate_testX.py` scripts predate the backend split and still call the
> deprecated flat API (`off.set_config`, `off.gen_nonbonded_tabpot`, …), which now emits
> `DeprecationWarning`. They still produce correct output. New tests should use
> `off.gmx.*` / `off.openmm.*`.

## Available Test Cases

### Test 1: Methane (Basic)
**Location**: `test/baseline_outputs/test1_methane_basic/`

**Purpose**: Tests basic workflow with no special options

**Covers**:
- Basic ReadOFF instantiation
- Nonbonded tabulated potential generation
- Bonded topology generation (minimal bonding)
- Default parameters

**Run**: `python generate_test1.py`

---

### Test 2: Ethane (Bonded Interactions)
**Location**: `test/baseline_outputs/test2_ethane_bonded/`

**Purpose**: Tests bonded interaction handling

**Covers**:
- Bond, angle, and dihedral interactions
- Bonded tabulated potential generation
- More complex molecular structure (C-C bond)

**Run**: `python generate_test2.py`

---

### Test 3: Water (Manual Charges)
**Location**: `test/baseline_outputs/test3_water_charges/`

**Purpose**: Tests that manually assigned charges appear correctly in topology

**Covers**:
- Manual charge assignment via `off.charges` dictionary
- Charge propagation to topology files
- `incl_mol` parameter usage
- Multiple molecule types in .off file (only H20QM used)

**Charges Used**: TIP3P-like (O=-0.82, H=+0.41, E=0.0)

**Run**: `python generate_test3.py`

**Note**: Charge calculation functions (calc_charges, normalization) have been removed from the project. This test demonstrates manual charge assignment, which is the only supported method.

---

### Test 4: Butanediol (Name Translation)
**Location**: `test/baseline_outputs/test4_butanediol_nametrans/`

**Purpose**: Tests atom name mapping between .off file and topology file

**Covers**:
- `name_translation` dictionary functionality
- Name mapping in tabulated potentials
- Name mapping in topology files

**Translation Map**:
```python
{
    'C1': 'CA',
    'C2': 'CB',
    'O1': 'OH',
    'H1': 'HO',
    'H2': 'HA',
    'H3': 'HB'
}
```

**Run**: `python generate_test4.py`

---

### Test 5: Methane (Soft-Core Sigma)
**Location**: `test/baseline_outputs/test5_methane_scsigma/`

**Purpose**: Tests soft-core scaling for free energy calculations

**Covers**:
- `sc_sigma` parameter (set to 0.3)
- Soft-core transformations in tabulated potentials
- Proper scaling for GROMACS free energy simulations

**Run**: `python generate_test5.py`

---

### Test 6: Ethane (Excluded Pairs)
**Location**: `test/baseline_outputs/test6_ethane_exclpairs/`

**Purpose**: Tests excluding specific atom pair interactions

**Covers**:
- `excl_pairs` parameter
- Pair exclusion in tabulated potentials
- Pair exclusion in topology files

**Excluded Pairs**: `[['C', 'C'], ['H', 'H']]` (only C-H interactions remain)

**Run**: `python generate_test6.py`

---

### Test 7: Methane (Configuration System)
**Location**: `test/baseline_outputs/test7_methane_config/`

**Purpose**: Tests configuration system features (set_config, get_config, load_charges_from_file)

**Covers**:
- `set_config()` method with multiple parameters
- `get_config()` method (all config and specific keys)
- `load_charges_from_file()` method
- Parameter resolution order (explicit → config → default)
- Method chaining for set_config() and load_charges_from_file()
- Explicit parameter override of config values

**Features Tested**:
- Configuration system correctly stores and retrieves values
- Workflow methods use config values when no explicit parameters provided
- Explicit parameters correctly override config values
- Charge file format with MOLNAME sections, comments, and blank lines
- Method chaining returns self for fluent API

**Run**: `python generate_test7.py`

---

### Test 8: Ethane (Clean Configuration-Based Workflow)
**Location**: `test/baseline_outputs/test8_ethane_clean_workflow/`

**Purpose**: Demonstrates clean, production-ready workflow using configuration system

**Covers**:
- Complete workflow from .off file to topology with config-based approach
- load_charges_from_file() in production context
- Bonded and nonbonded potential generation with config defaults
- Clean code pattern: set defaults once, use throughout

**Demonstrates**:
- How configuration system reduces code repetition
- Cleaner, more maintainable code compared to traditional approach
- Setting all defaults in one place
- Method calls without repetitive parameter passing

**Comparison**: Compare this test to Test 2 to see the difference between traditional parameter passing and config-based approach.

**Run**: `python generate_test8.py`

---

### Test 9: Butanol + Water (Reference FF Replacement)
**Location**: `test/test_regression.py::test9_change_molecule_blypsp4f` — no baseline directory

**Purpose**: Tests `change_molecule()` against the bundled `BLYPSP-4F` reference force field

**Covers**:
- Fitted water-water params (`EXPW`, `STRC`) are removed
- BLYPSP-4F water-water params (`EXP`, `STR`, `POW`) are inserted with correct values
- Solute-water **cross-term pairs keep their original fitted values**
- `H2OQM` charges are picked up from the sibling `BLYPSP-4F.charges`
- Solute (`UNK`) charges are unaffected

**Input**: `test/sample_off_files/h_butanol_fitwater.off`

**Run**: `python -m pytest test/test_regression.py::test9_change_molecule_blypsp4f -v`

This test is assertion-based rather than baseline-based — it checks parameter values
directly, so there is nothing to regenerate.

---

## How to Use the Regression Tests

### When Making Code Changes

1. **Before changing code**: confirm the suite is green
   ```bash
   python -m pytest test/ -v
   ```

2. **Make your code changes** (e.g., refactoring, adding features)

3. **Re-run the suite** — a failure means an output changed

4. **Interpret results**:
   - **All green**: ✅ your changes did not alter any output
   - **Failures**: investigate:
     - **Expected change**: you intentionally modified the output format → regenerate and
       commit new baselines (below)
     - **Unexpected change**: bug introduced → fix your code and re-test

Because `test_regression.py` writes into a `tmp_path`, running the suite never touches the
committed baselines — a failing test tells you something changed without overwriting the
evidence.

### Regenerating Baselines

The `generate_testN.py` scripts write directly into their own baseline directory, so run
them only when you *intend* to update baselines:

**Linux/Mac (Bash)**:
```bash
cd test/baseline_outputs
for dir in test*/; do
    (cd "$dir" && python generate_test*.py)
done
cd ../..
git diff test/baseline_outputs/       # review every change before committing
```

**Windows (PowerShell)**:
```powershell
cd test/baseline_outputs
foreach ($dir in Get-ChildItem -Directory) {
    cd $dir.Name
    python generate_test*.py
    cd ..
}
```

Then re-run `python -m pytest test/ -v` to confirm the new baselines are self-consistent.

### What Files to Check

After re-running tests, check these files for changes:
- `topol.top` - Final topology (most important)
- `temp_nonbonded.top` - Intermediate topology
- `tabpot/table*.xvg` - Tabulated potentials

**Note**: The generator scripts (`generate_test*.py`) should NOT change unless you're updating test logic.

## When to Update Baselines

Update baseline files (commit the changed outputs) when:

1. **Intentional output format change**: You modified how topology files are written
2. **Bug fix**: Previous baselines were incorrect due to a bug
3. **New feature with output impact**: Added feature that changes numerical output

**Never update baselines** when:
- Changes are unexpected (this indicates a bug)
- You don't understand why outputs changed
- Tests fail (fix the code first)

## Important Notes

### Template.top Format

All `template.top` files **must** have empty lines between bonded interaction sections:

```
[ atoms ]

[ bonds ]

[ angles ]

[ dihedrals ]

[ exclusions ]
```

Without empty lines, the topology generation will fail to populate these sections correctly.

### Charge Handling

**Important**: Charge calculation functions (calc_charges, charge normalization, residue neutralization) have been removed from the project. Only manual charge assignment is supported:

```python
off = afm.ReadOFF(off_loc="intra.off")
off.charges['MOLNAME']['ATOM1'] = -0.82
off.charges['MOLNAME']['ATOM2'] = 0.41
```

Test 3 demonstrates the correct approach for handling charges.

### Sample .off Files

Original sample .off files are located in `test/sample_off_files/`:
- `methane_intra.off` - Simple single carbon molecule
- `ethane_intra.off` - Small hydrocarbon with C-C bond
- `water_intra.off` - Multiple water models (H20QM, X3MM, etc.)
- `butanediol_intra.off` - 1,4-Butanediol (medium complexity)
- `big_alanine.off` - Large biomolecule (parsing tests only)
- `curcubituril.off` - Macrocycle (parsing tests only)
- `h_butanol_fitwater.off` - Butanol with co-fitted water (used by Test 9)

**Do not modify these files** - they are reference inputs for testing.

## Test Coverage Summary

| Feature | Test Coverage |
|---------|--------------|
| Basic workflow | Test 1 |
| Bonded interactions | Test 2, Test 8 |
| Manual charges | Test 3 |
| Name translation | Test 4 |
| Soft-core scaling | Test 5 |
| Pair exclusion | Test 6 |
| Configuration system | Test 7 |
| set_config() / get_config() | Test 7 |
| load_charges_from_file() | Test 7, Test 8 |
| Parameter resolution | Test 7 |
| Method chaining | Test 7 |
| Config-based workflow | Test 8 |
| C6 scaling | Test 1-8 (default) |
| Multiple molecules | Test 3 (H20QM only) |
| incl_mol parameter | Test 3 |
| change_molecule() / reference FFs | Test 9 |
| .off parsing (bonded + nonbonded) | `test_parsing.py`, 6 sample files |

## Not Covered

These have **no automated test** — changes to them are unguarded:

- `off.openmm.gen_xml()` and `off.openmm.preprocess_pdb()` — no OpenMM test exists
- `populate_bonded()` — no test, and no `atoms.dat` / `bonds.dat` / `valid_dihedrals.dat` /
  `parameters.dat` fixtures anywhere in the repo
- `write_report()`
- `compare.py`
- `gen_residues()`

## Future Test Additions

Consider adding tests for:
- The gaps listed under "Not Covered" above — a `populate_bonded` fixture directory is the
  highest-value one, since that code path is entirely unexercised
- `special_pairs` dictionary
- `excl_interactions` parameter
- `scale_C12` parameter (for HFE calculations)
- Multiple molecule types in same simulation
- Edge cases (empty molecules, missing parameters)

## Troubleshooting

### Test fails with ModuleNotFoundError
**Problem**: `ModuleNotFoundError: No module named 'numpy'`

**Solution**: Install numpy: `pip install numpy`

---

### Test fails with "no MOLNAME found in bonded_tabpot"
**Problem**: `incl_mol` parameter inconsistency

**Solution**: When using `incl_mol` in `gen_bonded_tabpot()`, you must also pass it to `gen_bonded_topology()` with the `bonded_tabpot` parameter:
```python
bonded_tabpot = off.gmx.gen_bonded_tabpot(incl_mol=['UNK'])
off.gmx.gen_bonded_topology(template_file='temp.top', write_to='topol.top',
                            incl_mol=['UNK'], bonded_tabpot=bonded_tabpot)
```

---

### Topology sections not populated
**Problem**: `[ atoms ]`, `[ bonds ]` sections are empty in `topol.top`

**Solution**: Ensure `template.top` has empty lines between sections (see "Template.top Format" above)

---

### Git shows unexpected changes
**Problem**: Re-running tests shows file differences even though code didn't change

**Solution**:
1. Check if the test generator script was modified
2. Verify numpy version (different versions may produce slightly different numerical output)
3. Check for non-deterministic behavior in code (random seeds, timestamps, etc.)

## Continuous Integration (Future)

No CI runner is configured yet. The suite is already CI-ready — it is plain pytest, needs
no installation step, and writes only into `tmp_path`:

```yaml
- run: pip install numpy pytest
- run: python -m pytest test/ -v
```

Adding `pytest-cov` (already in the `dev` extra) gives a coverage report:

```bash
python -m pytest test/ --cov=afmtogmx --cov-report=term-missing
```
