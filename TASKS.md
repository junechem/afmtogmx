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

**New signature** (parameter names match config keys, defaults to None):
```python
def gen_nonbonded_tabpot(self, special_pairs=None, incl_mol=None, excl_interactions=None,
                         excl_pairs=None, spacing_nonbonded=None, length_nonbonded=None,
                         scale_C6=None, sc_sigma=None)
```

**Add at the beginning of method** (before line 199):
```python
from types import SimpleNamespace

# Resolve parameters: explicit value → config → default (one line!)
p = SimpleNamespace(**{
    k: v if v is not None else self.config[k]
    for k, v in locals().items() if k in self.config
})

# Access via: p.spacing_nonbonded, p.incl_mol, p.special_pairs, etc.
# Update rest of function to use p.* instead of variable names
```

**Update docstring** to mention config fallback behavior and new parameter names.

---

### Task 3: Modify `gen_bonded_tabpot` Method
**File**: `src/afmtogmx/core/gen_md.py` (around line 217)

**Current signature**:
```python
def gen_bonded_tabpot(self, incl_mol=[], spacing=0.0001, length=0.3)
```

**New signature** (parameter names match config keys):
```python
def gen_bonded_tabpot(self, incl_mol=None, spacing_bonded=None, length_bonded=None)
```

**Add at the beginning of method**:
```python
from types import SimpleNamespace

# Resolve parameters: explicit value → config → default (one line!)
p = SimpleNamespace(**{
    k: v if v is not None else self.config[k]
    for k, v in locals().items() if k in self.config
})

# Access via: p.incl_mol, p.spacing_bonded, p.length_bonded
```

**Update docstring** to mention new parameter names.

---

### Task 4: Modify `gen_nonbonded_topology` Method
**File**: `src/afmtogmx/core/gen_md.py` (around line 241)

**Current signature**:
```python
def gen_nonbonded_topology(self, name_translation={}, template_file="", incl_mol=[],
                          excl_interactions=[], excl_pairs=[], write_to="", scale_C6=True,
                          scale_C12=1.0, special_pairs={}, sc_sigma=0.0)
```

**New signature** (parameter names match config keys):
```python
def gen_nonbonded_topology(self, name_translation=None, template_file="", incl_mol=None,
                          excl_interactions=None, excl_pairs=None, write_to="", scale_C6=None,
                          scale_C12=None, special_pairs=None, sc_sigma=None)
```

**Add at the beginning of method**:
```python
from types import SimpleNamespace

# Resolve parameters: explicit value → config → default (one line!)
p = SimpleNamespace(**{
    k: v if v is not None else self.config[k]
    for k, v in locals().items() if k in self.config
})

# Access via: p.name_translation, p.incl_mol, p.scale_C6, etc.
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

**New signature** (parameter names match config keys):
```python
def gen_bonded_topology(self, name_translation=None, incl_mol=None, template_file="",
                       write_to="", bonded_tabpot={})
```

**Add at the beginning of method**:
```python
from types import SimpleNamespace

# Resolve parameters: explicit value → config → default (one line!)
p = SimpleNamespace(**{
    k: v if v is not None else self.config[k]
    for k, v in locals().items() if k in self.config
})

# Access via: p.name_translation, p.incl_mol
```

**Update docstring**.

---

### Task 6: Modify `write_nonbonded_tabpot` Static Method
**File**: `src/afmtogmx/core/gen_md.py` (around line 359)

**Note**: This is currently a static method, which cannot access `self.config`.

**Approach**: Convert to instance method and keep static version for backward compatibility.

**Current signature**:
```python
@staticmethod
def write_nonbonded_tabpot(nonbonded_tabpot={}, prefix="MOL", to_dir="",
                          name_translation={}, write_blank=True):
```

**New implementation**:
```python
def write_nonbonded_tabpot(self, nonbonded_tabpot={}, tabpot_prefix=None, tabpot_dir=None,
                          name_translation=None, write_blank=None):
    """Instance method that uses config defaults."""
    from types import SimpleNamespace

    # Resolve parameters: explicit value → config → default (one line!)
    p = SimpleNamespace(**{
        k: v if v is not None else self.config[k]
        for k, v in locals().items() if k in self.config
    })

    # Call the static method with resolved parameters
    return ReadOFF._write_nonbonded_tabpot_static(nonbonded_tabpot, p.tabpot_prefix,
                                                   p.tabpot_dir, p.name_translation, p.write_blank)

@staticmethod
def _write_nonbonded_tabpot_static(nonbonded_tabpot={}, prefix="MOL", to_dir="",
                                   name_translation={}, write_blank=True):
    """Static implementation (for backward compatibility if called directly)."""
    # ... move existing implementation here ...
```

**Update docstring** to mention parameter name changes (`prefix` → `tabpot_prefix`, `to_dir` → `tabpot_dir`).

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
def write_bonded_tabpot(self, bonded_tabpot={}, tabpot_prefix=None, tabpot_dir=None):
    """Instance method that uses config defaults."""
    from types import SimpleNamespace

    # Resolve parameters: explicit value → config → default (one line!)
    p = SimpleNamespace(**{
        k: v if v is not None else self.config[k]
        for k, v in locals().items() if k in self.config
    })

    return ReadOFF._write_bonded_tabpot_static(bonded_tabpot, p.tabpot_prefix, p.tabpot_dir)

@staticmethod
def _write_bonded_tabpot_static(bonded_tabpot={}, prefix="MOL", to_dir=""):
    """Static implementation (for backward compatibility)."""
    # ... move existing implementation here ...
```

**Update docstring** to mention parameter name changes (`prefix` → `tabpot_prefix`, `to_dir` → `tabpot_dir`).

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
    tabpot = off.gen_nonbonded_tabpot(spacing_nonbonded=0.002)
    # Verify spacing_nonbonded=0.002 was used (check x-values in output)

def test_new_parameter_names():
    """Test new parameter names work correctly"""
    off = afm.ReadOFF(off_loc="test/sample_off_files/methane_intra.off")

    # New parameter names should work
    tabpot = off.gen_nonbonded_tabpot(
        excl_pairs=[['C1', 'H1']],
        spacing_nonbonded=0.001,
        length_nonbonded=3
    )
    # Verify it works

def test_method_chaining():
    """Test set_config returns self for chaining"""
    off = afm.ReadOFF(off_loc="test/sample_off_files/methane_intra.off")

    result = off.set_config({'sc_sigma': 0.5})
    assert result is off  # should return self
```

---

## Backward Compatibility Notes

**BREAKING CHANGES** (parameter name changes):
- `gen_nonbonded_tabpot`: `spacing` → `spacing_nonbonded`, `length` → `length_nonbonded`
- `gen_bonded_tabpot`: `spacing` → `spacing_bonded`, `length` → `length_bonded`
- `write_nonbonded_tabpot`: `prefix` → `tabpot_prefix`, `to_dir` → `tabpot_dir`
- `write_bonded_tabpot`: `prefix` → `tabpot_prefix`, `to_dir` → `tabpot_dir`

**Compatibility maintained**:
- [ ] Default behavior unchanged when parameters not provided (uses config defaults)
- [ ] Static method compatibility maintained (via `_static` wrappers)
- [ ] Empty/mutable default arguments handled correctly (None vs [] vs {})
- [ ] All tests must be updated to use new parameter names

---

## Migration Examples

### Before (Old Parameter Names - NO LONGER WORKS):
```python
off = afm.ReadOFF(off_loc="intra.off")
excl_pairs = [['OW', 'OW'], ['HW', 'EW']]
nonbonded_tabpot = off.gen_nonbonded_tabpot(excl_pairs=excl_pairs, sc_sigma=0.265, spacing=0.0005)
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, prefix='table')
off.gen_nonbonded_topology(template_file='template.top', write_to='temp.top',
                           excl_pairs=excl_pairs, sc_sigma=0.265)
off.gen_bonded_topology(template_file='temp.top', write_to='topol.top',
                       incl_mol=['UNK'])
```

### After (New Parameter Names - Without Config):
```python
off = afm.ReadOFF(off_loc="intra.off")
excl_pairs = [['OW', 'OW'], ['HW', 'EW']]
nonbonded_tabpot = off.gen_nonbonded_tabpot(excl_pairs=excl_pairs, sc_sigma=0.265, spacing_nonbonded=0.0005)
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, tabpot_prefix='table')
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

### Parameter Resolution Pattern

**Key Design Decision**: Parameter names match config keys exactly to enable automatic resolution.

**Implementation** (used in all methods):
```python
from types import SimpleNamespace

def method(self, spacing_nonbonded=None, incl_mol=None, ...):
    # One-line resolution: explicit value → config → default
    p = SimpleNamespace(**{
        k: v if v is not None else self.config[k]
        for k, v in locals().items() if k in self.config
    })

    # Use p.spacing_nonbonded, p.incl_mol, etc. throughout method
```

**Why this works**:
- `locals()` captures all function parameters and their values
- Filter by `k in self.config` to only process parameters that have config entries
- `SimpleNamespace` provides clean dot-notation access (`p.spacing_nonbonded`)
- No manual parameter listing needed - automatically picks up all config-aware parameters

### Handling Mutable Default Arguments

**IMPORTANT**: Using `None` as default avoids Python's mutable default argument issue where `[]` and `{}` are shared across calls.

### Config Precedence

Priority order (highest to lowest):
1. Explicitly provided parameter to method call (highest priority)
2. Value from `self.config` dictionary
3. Original hard-coded default value in `self.config`

### Method Chaining

`set_config` returns `self` to enable fluent API:
```python
off = afm.ReadOFF("intra.off").set_config({...}).load_charges_from_file(...)
```

### Parameter Naming Convention

Parameters that map to config use the **exact config key name**:
- `spacing_nonbonded` (not `spacing`)
- `length_nonbonded` (not `length`)
- `tabpot_prefix` (not `prefix`)
- `tabpot_dir` (not `to_dir`)

This enables automatic parameter resolution without manual mapping.

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
