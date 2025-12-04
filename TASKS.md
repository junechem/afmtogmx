# TASKS.md

## Completed Tasks

### ✅ Project Setup & Documentation (2025-12-03)
- [x] Set up GitHub repository (junechem/afmtogmx)
- [x] Configure git credentials for junechem account
- [x] Create CLAUDE.md with comprehensive project documentation
- [x] Create environment.yml for conda dependencies
- [x] Create pyproject.toml for Python package configuration
- [x] Update .gitignore for Python projects
- [x] Create main and in_progress branches

### ✅ Regression Testing Infrastructure (2025-12-04)
- [x] Create test/baseline_outputs/ directory structure
- [x] Install numpy dependency
- [x] Create Test 1: Methane (basic) - tests basic workflow with no special options
- [x] Create Test 2: Ethane (bonded) - tests bonded interactions (bonds, angles, dihedrals)
- [x] Create Test 3: Water (manual charges) - tests manually assigned charges in topology
- [x] Create Test 4: Butanediol (name_translation) - tests atom name mapping
- [x] Create Test 5: Methane (sc_sigma) - tests soft-core scaling for free energy calculations
- [x] Create Test 6: Ethane (excl_pairs) - tests excluding specific atom pair interactions
- [x] Document testing approach in TESTING.md
- **Note**: All template.top files need empty lines between [ atoms ], [ bonds ], [ angles ], etc. sections

### ✅ Remove Charge Calculation Functions (2025-12-04)
- [x] Removed calc_charges() method from ReadOFF class
- [x] Deleted chargefxns.py module entirely (192 lines)
- [x] Removed charge-only functions from residues.py (_check_molname_resname, _generate_residue_priority)
- [x] Updated __init__.py documentation to show manual charge assignment
- [x] Updated CLAUDE.md to remove all charge calculation references
- [x] Verified all 6 regression tests still pass
- **Result**: Only manual charge assignment via off.charges dictionary is now supported
- **Net change**: Removed 286 lines of code (314 deletions, 28 insertions)

## In Progress

(No active tasks)

## Future Work

### Configuration Object Pattern Implementation

#### Overview
Implement a configuration system for the ReadOFF class to avoid repeating common parameters across multiple method calls. This maintains backward compatibility while providing a cleaner API.

#### Goals
- Store shared configuration options once in the ReadOFF instance
- Allow methods to use config as fallback when parameters not explicitly provided
- Maintain 100% backward compatibility with existing code
- Enable method chaining for fluent API

## Implementation Tasks

### Task 1: Add Configuration Infrastructure to ReadOFF Class
**File**: `src/afmtogmx/core/gen_md.py`

**Changes to `__init__` method** (around line 10-40):
- Add `self.config = {}` dictionary after existing instance variables
- Initialize with default configuration values:
  ```python
  self.config = {
      'incl_mol': [],
      'excl_pairs': [],
      'excl_interactions': [],
      'special_pairs': {},
      'name_translation': {},
      'sc_sigma': 0.0,
      'spacing_nonbonded': 0.0005,
      'spacing_bonded': 0.0001,
      'length_nonbonded': 3,
      'length_bonded': 0.3,
      'scale_C6': True,
      'scale_C12': 1.0,
      'tabpot_prefix': 'MOL',
      'tabpot_dir': '',
      'write_blank': True,
  }
  ```

**Add `set_config` method** (after `__init__`):
```python
def set_config(self, config_dict):
    """Update configuration options that will be used as defaults across methods.

    This allows setting common parameters once instead of repeating them in every
    method call. Individual method calls can still override config values.

    :param config_dict: Dictionary with configuration keys and values
    :return: self (for method chaining)

    Example:
        >>> off = ReadOFF(off_loc="intra.off")
        >>> off.set_config({
        ...     'incl_mol': ['UNK'],
        ...     'excl_pairs': [['OW', 'OW'], ['HW', 'EW']],
        ...     'sc_sigma': 0.265,
        ... })
        >>> nonbonded_tabpot = off.gen_nonbonded_tabpot()  # uses config values
    """
    self.config.update(config_dict)
    return self
```

**Add `get_config` method** (helper for debugging):
```python
def get_config(self, key=None):
    """Get current configuration value(s).

    :param key: Optional specific config key. If None, returns entire config dict
    :return: config value or entire config dict
    """
    if key is None:
        return self.config.copy()
    return self.config.get(key)
```

---

### Task 2: Modify `gen_nonbonded_tabpot` Method
**File**: `src/afmtogmx/core/gen_md.py` (around line 137)

**Current signature**:
```python
def gen_nonbonded_tabpot(self, special_pairs={}, incl_mol=[], excl_interactions=[],
                         excl_pairs=[], spacing=0.0005, length=3, scale_C6=True, sc_sigma=0.0)
```

**New signature** (change defaults to None):
```python
def gen_nonbonded_tabpot(self, special_pairs=None, incl_mol=None, excl_interactions=None,
                         excl_pairs=None, spacing=None, length=None, scale_C6=None, sc_sigma=None)
```

**Add at the beginning of method** (before line 199):
```python
# Use provided values, or fall back to config, or use original defaults
special_pairs = special_pairs if special_pairs is not None else self.config.get('special_pairs', {})
incl_mol = incl_mol if incl_mol is not None else self.config.get('incl_mol', [])
excl_interactions = excl_interactions if excl_interactions is not None else self.config.get('excl_interactions', [])
excl_pairs = excl_pairs if excl_pairs is not None else self.config.get('excl_pairs', [])
spacing = spacing if spacing is not None else self.config.get('spacing_nonbonded', 0.0005)
length = length if length is not None else self.config.get('length_nonbonded', 3)
scale_C6 = scale_C6 if scale_C6 is not None else self.config.get('scale_C6', True)
sc_sigma = sc_sigma if sc_sigma is not None else self.config.get('sc_sigma', 0.0)
```

**Update docstring** to mention config fallback behavior.

---

### Task 3: Modify `gen_bonded_tabpot` Method
**File**: `src/afmtogmx/core/gen_md.py` (around line 217)

**Current signature**:
```python
def gen_bonded_tabpot(self, incl_mol=[], spacing=0.0001, length=0.3)
```

**New signature**:
```python
def gen_bonded_tabpot(self, incl_mol=None, spacing=None, length=None)
```

**Add at the beginning of method**:
```python
# Use provided values, or fall back to config, or use original defaults
incl_mol = incl_mol if incl_mol is not None else self.config.get('incl_mol', [])
spacing = spacing if spacing is not None else self.config.get('spacing_bonded', 0.0001)
length = length if length is not None else self.config.get('length_bonded', 0.3)
```

**Update docstring**.

---

### Task 4: Modify `gen_nonbonded_topology` Method
**File**: `src/afmtogmx/core/gen_md.py` (around line 241)

**Current signature**:
```python
def gen_nonbonded_topology(self, name_translation={}, template_file="", incl_mol=[],
                          excl_interactions=[], excl_pairs=[], write_to="", scale_C6=True,
                          scale_C12=1.0, special_pairs={}, sc_sigma=0.0)
```

**New signature**:
```python
def gen_nonbonded_topology(self, name_translation=None, template_file="", incl_mol=None,
                          excl_interactions=None, excl_pairs=None, write_to="", scale_C6=None,
                          scale_C12=None, special_pairs=None, sc_sigma=None)
```

**Add at the beginning of method**:
```python
# Use provided values, or fall back to config, or use original defaults
name_translation = name_translation if name_translation is not None else self.config.get('name_translation', {})
incl_mol = incl_mol if incl_mol is not None else self.config.get('incl_mol', [])
excl_interactions = excl_interactions if excl_interactions is not None else self.config.get('excl_interactions', [])
excl_pairs = excl_pairs if excl_pairs is not None else self.config.get('excl_pairs', [])
scale_C6 = scale_C6 if scale_C6 is not None else self.config.get('scale_C6', True)
scale_C12 = scale_C12 if scale_C12 is not None else self.config.get('scale_C12', 1.0)
special_pairs = special_pairs if special_pairs is not None else self.config.get('special_pairs', {})
sc_sigma = sc_sigma if sc_sigma is not None else self.config.get('sc_sigma', 0.0)
```

**Update docstring**.

---

### Task 5: Modify `gen_bonded_topology` Method
**File**: `src/afmtogmx/core/gen_md.py` (around line 304)

**Current signature**:
```python
def gen_bonded_topology(self, name_translation={}, incl_mol=[], template_file="",
                       write_to="", bonded_tabpot={})
```

**New signature**:
```python
def gen_bonded_topology(self, name_translation=None, incl_mol=None, template_file="",
                       write_to="", bonded_tabpot={})
```

**Add at the beginning of method**:
```python
# Use provided values, or fall back to config, or use original defaults
name_translation = name_translation if name_translation is not None else self.config.get('name_translation', {})
incl_mol = incl_mol if incl_mol is not None else self.config.get('incl_mol', [])
```

**Update docstring**.

---

### Task 6: Modify `write_nonbonded_tabpot` Static Method
**File**: `src/afmtogmx/core/gen_md.py` (around line 359)

**Note**: This is currently a static method, which cannot access `self.config`.

**Options**:
1. Make it an instance method (recommended)
2. Keep it static but add a convenience instance method wrapper

**Recommended approach**: Convert to instance method and add wrapper for backward compatibility.

**Current signature**:
```python
@staticmethod
def write_nonbonded_tabpot(nonbonded_tabpot={}, prefix="MOL", to_dir="",
                          name_translation={}, write_blank=True):
```

**New implementation**:
```python
def write_nonbonded_tabpot(self, nonbonded_tabpot={}, prefix=None, to_dir=None,
                          name_translation=None, write_blank=None):
    """Instance method that uses config defaults."""
    # Use provided values, or fall back to config, or use original defaults
    prefix = prefix if prefix is not None else self.config.get('tabpot_prefix', 'MOL')
    to_dir = to_dir if to_dir is not None else self.config.get('tabpot_dir', '')
    name_translation = name_translation if name_translation is not None else self.config.get('name_translation', {})
    write_blank = write_blank if write_blank is not None else self.config.get('write_blank', True)

    # Call the static method with resolved parameters
    return ReadOFF._write_nonbonded_tabpot_static(nonbonded_tabpot, prefix, to_dir,
                                                   name_translation, write_blank)

@staticmethod
def _write_nonbonded_tabpot_static(nonbonded_tabpot={}, prefix="MOL", to_dir="",
                                   name_translation={}, write_blank=True):
    """Static implementation (for backward compatibility if called directly)."""
    # ... move existing implementation here ...
```

**Update docstring**.

---

### Task 7: Modify `write_bonded_tabpot` Static Method
**File**: `src/afmtogmx/core/gen_md.py` (around line 415)

**Current signature**:
```python
@staticmethod
def write_bonded_tabpot(bonded_tabpot={}, prefix="MOL", to_dir=""):
```

**New implementation** (same pattern as Task 6):
```python
def write_bonded_tabpot(self, bonded_tabpot={}, prefix=None, to_dir=None):
    """Instance method that uses config defaults."""
    prefix = prefix if prefix is not None else self.config.get('tabpot_prefix', 'MOL')
    to_dir = to_dir if to_dir is not None else self.config.get('tabpot_dir', '')
    return ReadOFF._write_bonded_tabpot_static(bonded_tabpot, prefix, to_dir)

@staticmethod
def _write_bonded_tabpot_static(bonded_tabpot={}, prefix="MOL", to_dir=""):
    """Static implementation (for backward compatibility)."""
    # ... move existing implementation here ...
```

**Update docstring**.

---

### Task 8: Add Convenience Helper - Load Charges from File
**File**: `src/afmtogmx/core/gen_md.py`

**Add new method** (after `calc_charges`):
```python
def load_charges_from_file(self, filepath, molname, delimiter=None, comment='#'):
    """Load atomic charges from a text file.

    :param filepath: Path to file containing charges
    :param molname: Molecule name to assign charges to
    :param delimiter: Column delimiter (default: any whitespace)
    :param comment: Comment character to ignore (default: '#')

    File format (two columns: atom_name charge):
        OW  -0.82
        HW   0.41
        # This is a comment
        EW   0.00
    """
    if molname not in self.charges:
        print(f"Warning: molname '{molname}' not found in charges dictionary")
        return self

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(comment):
                continue

            parts = line.split(delimiter)
            if len(parts) < 2:
                continue

            atom = parts[0]
            try:
                charge = float(parts[1])
                if atom in self.charges[molname]:
                    self.charges[molname][atom] = charge
                else:
                    print(f"Warning: atom '{atom}' not found in molecule '{molname}'")
            except ValueError:
                print(f"Warning: could not parse charge from line: {line}")

    return self  # for chaining
```

---

### Task 9: Update Package-Level Documentation
**File**: `src/afmtogmx/__init__.py`

**Add new example workflow** (after existing examples, around line 68):
```python
Alternatively, you can use the configuration system to avoid repeating parameters:

    >>> import afmtogmx as afm
    >>> off = afm.ReadOFF(off_loc="path/to/intra.off")
    >>> off.set_config({
    ...     'incl_mol': ['UNK'],
    ...     'excl_pairs': [['OW', 'OW'], ['HW', 'EW']],
    ...     'sc_sigma': 0.265,
    ...     'tabpot_prefix': 'table'
    ... })
    >>> off.load_charges_from_file('final-charges', molname='UNK')
    >>> nonbonded_tabpot = off.gen_nonbonded_tabpot()
    >>> off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot)
    >>> off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')
    >>> off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')
```

---

### Task 10: Update CLAUDE.md Documentation
**File**: `CLAUDE.md`

**Add new section** after "Standard Workflow" (around line 95):
```markdown
## Configuration System (Recommended)

To avoid repeating common parameters, use the configuration system:

```python
off = afm.ReadOFF(off_loc="path/to/intra.off")

# Set configuration once
off.set_config({
    'incl_mol': ['UNK'],
    'excl_pairs': [['OW', 'OW'], ['HW', 'EW']],
    'sc_sigma': 0.265,
    'tabpot_prefix': 'table',
})

# Methods use config values as defaults (but can still override)
nonbonded_tabpot = off.gen_nonbonded_tabpot()  # uses config
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot)

# Override config for specific call
bonded_tabpot = off.gen_bonded_tabpot(spacing=0.00005)  # override spacing
```

### Configuration Keys

Available configuration options:
- `incl_mol`: List of molecule names to include (default: [])
- `excl_pairs`: List of atom pairs to exclude (default: [])
- `excl_interactions`: List of interaction names to exclude (default: [])
- `special_pairs`: Dict of special pair treatments (default: {})
- `name_translation`: Dict for atom name mapping (default: {})
- `sc_sigma`: Soft-core sigma for FE calculations (default: 0.0)
- `spacing_nonbonded`: Spacing for nonbonded tables (default: 0.0005)
- `spacing_bonded`: Spacing for bonded tables (default: 0.0001)
- `length_nonbonded`: Length for nonbonded tables (default: 3)
- `length_bonded`: Length for bonded tables (default: 0.3)
- `scale_C6`: Scale C6 for dispersion (default: True)
- `scale_C12`: Scale C12 for HFE (default: 1.0)
- `tabpot_prefix`: Prefix for tabulated potential files (default: 'MOL')
- `tabpot_dir`: Directory for tabulated potentials (default: '')
- `write_blank`: Write blank tables (default: True)
```

**Update "Recent Changes" section** to mention the configuration system addition.

---

## Testing Requirements

### Test File: Create `test/test_config.py`

**Test cases needed**:
1. Test `set_config` updates config dictionary
2. Test methods use config values when parameters not provided
3. Test explicit parameters override config values
4. Test backward compatibility (old code still works)
5. Test method chaining: `off.set_config({...}).gen_nonbonded_tabpot()`
6. Test `load_charges_from_file` helper
7. Test `get_config` retrieves values correctly
8. Test config doesn't affect different ReadOFF instances

**Example test structure**:
```python
import afmtogmx as afm

def test_config_basic():
    """Test basic config set/get"""
    off = afm.ReadOFF(off_loc="test/sample_off_files/methane_intra.off")
    off.set_config({'incl_mol': ['TEST'], 'sc_sigma': 0.5})

    assert off.get_config('incl_mol') == ['TEST']
    assert off.get_config('sc_sigma') == 0.5

def test_config_fallback():
    """Test methods use config as fallback"""
    off = afm.ReadOFF(off_loc="test/sample_off_files/methane_intra.off")
    off.set_config({'excl_pairs': [['C1', 'H1']], 'spacing_nonbonded': 0.001})

    # Would need to inspect that methods actually use these values
    # This requires either returning the used config or checking outputs

def test_config_override():
    """Test explicit params override config"""
    off = afm.ReadOFF(off_loc="test/sample_off_files/methane_intra.off")
    off.set_config({'spacing_nonbonded': 0.001})

    # Call with explicit spacing should override config
    tabpot = off.gen_nonbonded_tabpot(spacing=0.002)
    # Verify spacing=0.002 was used (check x-values in output)

def test_backward_compatibility():
    """Test old code still works without config"""
    off = afm.ReadOFF(off_loc="test/sample_off_files/methane_intra.off")

    # Old-style calls should work exactly as before
    tabpot = off.gen_nonbonded_tabpot(
        excl_pairs=[['C1', 'H1']],
        spacing=0.001
    )
    # Verify it works

def test_method_chaining():
    """Test set_config returns self for chaining"""
    off = afm.ReadOFF(off_loc="test/sample_off_files/methane_intra.off")

    result = off.set_config({'sc_sigma': 0.5})
    assert result is off  # should return self
```

---

## Backward Compatibility Checklist

- [ ] All existing code continues to work without modifications
- [ ] Default behavior unchanged when config not used
- [ ] Static method compatibility maintained (via _static wrappers)
- [ ] No breaking changes to method signatures (only default values changed)
- [ ] Empty/mutable default arguments handled correctly (None vs [] vs {})

---

## Migration Examples

### Before (Current API):
```python
off = afm.ReadOFF(off_loc="intra.off")
excl_pairs = [['OW', 'OW'], ['HW', 'EW']]
nonbonded_tabpot = off.gen_nonbonded_tabpot(excl_pairs=excl_pairs, sc_sigma=0.265)
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, prefix='table')
off.gen_nonbonded_topology(template_file='template.top', write_to='temp.top',
                           excl_pairs=excl_pairs, sc_sigma=0.265)
off.gen_bonded_topology(template_file='temp.top', write_to='topol.top',
                       incl_mol=['UNK'])
```

### After (With Config):
```python
off = afm.ReadOFF(off_loc="intra.off")
off.set_config({
    'incl_mol': ['UNK'],
    'excl_pairs': [['OW', 'OW'], ['HW', 'EW']],
    'sc_sigma': 0.265,
    'tabpot_prefix': 'table'
})
nonbonded_tabpot = off.gen_nonbonded_tabpot()
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot)
off.gen_nonbonded_topology(template_file='template.top', write_to='temp.top')
off.gen_bonded_topology(template_file='temp.top', write_to='topol.top')
```

---

## Implementation Notes

### Handling Mutable Default Arguments

**IMPORTANT**: Python's mutable default arguments (`[]`, `{}`) are shared across calls. Using `None` as default and replacing with actual defaults in the method body avoids this issue.

**Bad**:
```python
def method(self, excl_pairs=[]):  # Same list reused!
    excl_pairs.append(something)
```

**Good**:
```python
def method(self, excl_pairs=None):
    excl_pairs = excl_pairs if excl_pairs is not None else self.config.get('excl_pairs', [])
    excl_pairs = excl_pairs.copy()  # or list(excl_pairs) to avoid modifying config
```

### Config Precedence

Priority order (highest to lowest):
1. Explicitly provided parameter to method call
2. Value from `self.config` dictionary
3. Original hard-coded default value

### Method Chaining

`set_config` returns `self` to enable:
```python
off = afm.ReadOFF("intra.off").set_config({...}).load_charges_from_file(...)
```

---

## Estimated Effort

- **Task 1**: 15 minutes (add config infrastructure)
- **Task 2-7**: 10 minutes each (60 minutes total - modify methods)
- **Task 8**: 20 minutes (add load_charges helper)
- **Task 9-10**: 15 minutes (update documentation)
- **Testing**: 45 minutes (write and run tests)

**Total: ~2.5 hours**

---

## Success Criteria

- [ ] All methods accept config values as defaults
- [ ] `set_config()` and `get_config()` work correctly
- [ ] Method chaining works: `off.set_config({...}).gen_nonbonded_tabpot()`
- [ ] User example workflow reduced from repeating parameters 3-4 times to once
- [ ] All existing code continues to work (backward compatible)
- [ ] Tests pass
- [ ] Documentation updated (CLAUDE.md, __init__.py)
- [ ] `load_charges_from_file` helper implemented
