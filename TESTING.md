# Testing Guide

This document explains the regression testing infrastructure for the `afmtogmx` project.

## Overview

The regression testing system ensures that code changes don't break existing functionality. Tests are located in `test/baseline_outputs/` and consist of baseline output files that represent correct behavior.

## Test Structure

Each test directory contains:
- `generate_testX.py` - Python script to regenerate outputs
- `template.top` - Input GROMACS topology template
- `topol.top` - Final topology file (baseline)
- `temp_nonbonded.top` - Intermediate topology file (baseline)
- `tabpot/*.xvg` - Tabulated potential files (baseline)

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

**Note**: Charge calculation functions (calc_charges, normalization) will be removed from the project. This test demonstrates manual charge assignment, which is the only supported method going forward.

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

## How to Use Regression Tests

### When Making Code Changes

1. **Before changing code**: Ensure all current tests pass
   ```bash
   cd test/baseline_outputs/test1_methane_basic && python generate_test1.py
   cd ../test2_ethane_bonded && python generate_test2.py
   cd ../test3_water_charges && python generate_test3.py
   cd ../test4_butanediol_nametrans && python generate_test4.py
   cd ../test5_methane_scsigma && python generate_test5.py
   cd ../test6_ethane_exclpairs && python generate_test6.py
   ```

2. **Make your code changes** (e.g., refactoring, adding features)

3. **Re-run all test generators** (same commands as step 1)

4. **Check for differences**:
   ```bash
   git status  # Shows which files changed
   git diff test/baseline_outputs/  # Shows exact line-by-line differences
   ```

5. **Interpret results**:
   - **No changes**: ✅ Your code changes are safe and backward-compatible
   - **Changes detected**: Investigate further:
     - **Expected changes**: You intentionally modified output format → Commit new baselines
     - **Unexpected changes**: Bug introduced → Fix your code and re-test

### Running All Tests Quickly

You can create a simple script to run all tests:

**Windows (PowerShell)**:
```powershell
cd test/baseline_outputs
foreach ($dir in Get-ChildItem -Directory) {
    cd $dir.Name
    python generate_test*.py
    cd ..
}
```

**Linux/Mac (Bash)**:
```bash
cd test/baseline_outputs
for dir in test*/; do
    cd "$dir"
    python generate_test*.py
    cd ..
done
```

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

**Important**: Charge calculation functions (calc_charges, charge normalization, residue neutralization) are being removed from the project. Only manual charge assignment is supported:

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
- `big_alanine.off` - Large biomolecule (not used in current tests)

**Do not modify these files** - they are reference inputs for testing.

## Test Coverage Summary

| Feature | Test Coverage |
|---------|--------------|
| Basic workflow | Test 1 |
| Bonded interactions | Test 2 |
| Manual charges | Test 3 |
| Name translation | Test 4 |
| Soft-core scaling | Test 5 |
| Pair exclusion | Test 6 |
| C6 scaling | Test 1-6 (default) |
| Multiple molecules | Test 3 (H20QM only) |
| incl_mol parameter | Test 3 |

## Future Test Additions

Consider adding tests for:
- `special_pairs` dictionary
- `excl_interactions` parameter
- `scale_C12` parameter (for HFE calculations)
- Large molecules (big_alanine.off)
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
bonded_tabpot = off.gen_bonded_tabpot(incl_mol=['UNK'])
off.gen_bonded_topology(template_file='temp.top', write_to='topol.top',
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

For CI/CD integration, create a test runner that:
1. Runs all test generators
2. Checks git diff for changes
3. Fails if unexpected changes detected
4. Reports which tests passed/failed

Example structure:
```python
import subprocess
import sys

tests = ['test1_methane_basic', 'test2_ethane_bonded', ...]
failed = []

for test in tests:
    result = subprocess.run(['python', f'test/baseline_outputs/{test}/generate_{test}.py'])
    if result.returncode != 0:
        failed.append(test)

# Check for git changes
diff_result = subprocess.run(['git', 'diff', '--exit-code', 'test/baseline_outputs/'])
if diff_result.returncode != 0:
    print("ERROR: Baseline outputs changed unexpectedly!")
    sys.exit(1)
```
