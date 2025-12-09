# TASKS.md

## In Progress

(No active tasks)

## Completed

### ✅ Simplify Tabulated Potential Workflow with Attribute Storage (2025-12-09)

**Status**: COMPLETED - All tests passing ✓

**Goal**: Eliminate the need to store tabulated potentials in variables by storing them as attributes on the ReadOFF object.

**Before**:
```python
nonbonded_tabpot = off.gen_nonbonded_tabpot()  # Must save to variable
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot)  # Pass it back
```

**After**:
```python
off.gen_nonbonded_tabpot()  # Generates and stores internally
off.write_nonbonded_tabpot()  # Uses stored value automatically
```

**Implementation Summary**:

✅ All 9 tasks completed successfully:

1. ✅ Updated ReadOFF.__init__() - Added `self.nonbonded_tabpot` and `self.bonded_tabpot` attributes
2. ✅ Updated gen_nonbonded_tabpot() - Stores result in `self.nonbonded_tabpot`
3. ✅ Updated gen_bonded_tabpot() - Stores result in `self.bonded_tabpot`
4. ✅ Updated write_nonbonded_tabpot() - Uses `self.nonbonded_tabpot` by default
5. ✅ Updated write_bonded_tabpot() - Uses `self.bonded_tabpot` by default
6. ✅ Updated test8 to demonstrate new simplified workflow
7. ✅ Updated CLAUDE.md documentation with new workflow examples
8. ✅ Updated __init__.py workflow examples
9. ✅ All regression tests pass (tests 1-8)

**Files Modified**:
- `src/afmtogmx/core/gen_md.py`: Core implementation
- `src/afmtogmx/__init__.py`: Workflow documentation
- `CLAUDE.md`: Main documentation
- `test/baseline_outputs/test8_ethane_clean_workflow/generate_test8.py`: Simplified workflow example

**Key Features**:
- Maintains 100% backward compatibility (can still pass parameters explicitly)
- Still allows flexibility to generate multiple variants and pass custom dicts
- Attributes get overwritten each time gen methods are called (acceptable limitation)
- Methods still return dictionaries for use cases that need to capture/inspect them
- Clear error messages if write methods called before gen methods

**Test Results**:
```
✓ test1_methane_basic PASSED
✓ test2_ethane_bonded PASSED (backward compatibility verified)
✓ test3_water_charges PASSED
✓ test4_butanediol_nametrans PASSED
✓ test5_methane_scsigma PASSED
✓ test6_ethane_exclpairs PASSED
✓ test7_methane_config PASSED
✓ test8_ethane_clean_workflow PASSED (new simplified workflow)
```

## Future Work

(Add future tasks here)
