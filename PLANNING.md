# PLANNING.md — OpenMM Integration & Backend Architecture

## Motivation

`afmtogmx` (GROMACS outputs) and `offtoopenmm` (OpenMM outputs) both parse the same
CRYOFF `.off` files. Rather than fixing the 5 architectural problems in `offtoopenmm`
in isolation, the right solution is to merge OpenMM output into `afmtogmx` so both
engines share one parser.

The key practical benefit: when an AFM-fitted `.off` file contains refitted water-water
parameters (which you don't want), you can call `off.change_molecule('H2OQM', 'BLYPSP-4F')`
to replace them in-place with the stored reference FF before generating any output.
This produces a complete, self-contained OpenMM XML — no need to load multiple XML files
at runtime, no naming conflicts, no matrix-exclusion complexity.

---

## New API

### Shared on `ReadOFF` (unchanged behavior)

```python
off = afm.ReadOFF('forcefield.off')

# Shared data — populated at parse time
off.bonded
off.nonbonded
off.charges
off.residues

# Shared methods — modify the above data structures
off.load_charges_from_file('charges.txt')
off.gen_residues(residue_definition={}, residue_atnums={})

# NEW shared method — replaces a molecule's parameters with a stored reference FF
off.change_molecule(
    mol_name='H2OQM',               # molecule name as it appears in the .off file
    reference_ff='BLYPSP-4F',       # name of stored reference FF
    atom_name_map={'OW': 'OW_sp4f', 'HW': 'HW_sp4f'}  # rename atom types to match reference FF
)

# Engine-specific backends
off.gmx     # GROMACSBackend
off.openmm  # OpenMMBackend
```

### `off.gmx` — GROMACS backend

```python
off.gmx.set_config(
    spacing_nonbonded=0.0005,
    length_nonbonded=3.0,
    spacing_bonded=0.0001,
    length_bonded=0.3,
    tabpot_prefix='table',
    tabpot_dir='',
    scale_C6=True,
    sc_sigma=0.0,
    write_blank=True,
    name_translation={},
    special_pairs={},
    incl_mol=[],
    excl_interactions=[],
    excl_pairs=[],
)
off.gmx.get_config()               # returns full config dict
off.gmx.get_config('scale_C6')    # returns specific value

off.gmx.gen_nonbonded_tabpot(...)  # stores in off.gmx.nonbonded_tabpot
off.gmx.write_nonbonded_tabpot(...)
off.gmx.gen_bonded_tabpot(...)     # stores in off.gmx.bonded_tabpot
off.gmx.write_bonded_tabpot(...)
off.gmx.gen_nonbonded_topology(template_file='template.top', write_to='nonbonded.top')
off.gmx.gen_bonded_topology(template_file='nonbonded.top', write_to='topol.top')
```

### `off.openmm` — OpenMM backend

```python
off.openmm.set_config(
    incl_mol=[],
    molname_translations={},   # e.g. {'H2OQM': 'SOL', 'UNK': 'UNK'}
)
off.openmm.get_config()

off.openmm.gen_xml(output='forcefield.xml', ...)
```

### Deprecated — old API still works, emits `DeprecationWarning`

All old methods on `ReadOFF` become thin wrappers that warn and delegate to `off.gmx`:

```python
off.set_config(...)           # → off.gmx.set_config(...)
off.get_config(...)           # → off.gmx.get_config(...)
off.gen_nonbonded_tabpot(...) # → off.gmx.gen_nonbonded_tabpot(...)
off.write_nonbonded_tabpot(...)
off.gen_bonded_tabpot(...)
off.write_bonded_tabpot(...)
off.gen_nonbonded_topology(...)
off.gen_bonded_topology(...)

# These become properties that forward to off.gmx equivalents
off.config            # → off.gmx.config
off.nonbonded_tabpot  # → off.gmx.nonbonded_tabpot
off.bonded_tabpot     # → off.gmx.bonded_tabpot
```

Deprecated methods return `self` (the `ReadOFF` object, not `off.gmx`) so old
method chaining continues to work without changes.

---

## Package Structure

Minimal restructuring — add new files, don't move existing ones.

```
src/afmtogmx/
├── __init__.py
├── core/
│   ├── gen_md.py              — ReadOFF: parsing + shared methods + deprecation wrappers
│   ├── gmx_backend.py         — NEW: GROMACSBackend class
│   ├── openmm_backend.py      — NEW: OpenMMBackend class
│   ├── functions.py           — unchanged
│   ├── tabulated_potentials.py — unchanged (used by gmx_backend.py)
│   ├── topology.py            — unchanged (used by gmx_backend.py)
│   ├── residues.py            — unchanged
│   └── compare.py             — unchanged
└── reference_ff/
    └── BLYPSP-4F.off          — NEW: stored reference FF (partial .off format)
```

`pyproject.toml` packages list will need updating once `reference_ff/` is added.

---

## `change_molecule()` — How It Works

```python
off.change_molecule(
    mol_name='H2OQM',
    reference_ff='BLYPSP-4F',
    atom_name_map={'OW': 'OW_blys', 'HW': 'HW_blys'}
)
```

1. Load `src/afmtogmx/reference_ff/BLYPSP-4F.off` using the existing parser
2. Find all atom types belonging to `mol_name` in `off.nonbonded`, `off.bonded`, `off.charges`
3. Rename those atom types using `atom_name_map` (old name → new name matching reference FF)
4. Overwrite parameter values with those from the reference FF

The atom name mapping is necessary because the AFM-fitted `.off` file may use different
atom type names for water than the stored reference FF does.

---

## Implementation Order

1. ✅ **`gmx_backend.py`** — Create `GROMACSBackend` class. Move all GROMACS methods
   (`gen_nonbonded_tabpot`, `write_nonbonded_tabpot`, `gen_bonded_tabpot`,
   `write_bonded_tabpot`, `gen_nonbonded_topology`, `gen_bonded_topology`,
   `set_config`, `get_config`) from `gen_md.py` into this class.
   Each method reads from `self._parent.nonbonded`, `self._parent.bonded`, etc.

2. ✅ **`ReadOFF.__init__`** — Add `self.gmx = GROMACSBackend(self)`.
   Remove old method bodies (they become wrappers).

3. ✅ **Deprecation wrappers** — Add deprecated methods and properties to `ReadOFF`
   that warn and delegate to `off.gmx` (`config`, `nonbonded_tabpot`, `bonded_tabpot`,
   `set_config`, `get_config`, and all six generation/write methods).

4. ✅ **`openmm_backend.py`** — Create `OpenMMBackend` skeleton with `gen_xml()` stub
   (raises `NotImplementedError` until Step 8).

5. **Wire `OpenMMBackend`** — Add `self.openmm = OpenMMBackend(self)` to
   `ReadOFF.__init__`.

6. **Update tests** — Update `test_regression.py` to use `off.gmx.xxx()`.
   `test_parsing.py` needs no changes (tests shared data structures). Do this
   before the OpenMM port to confirm the GROMACS path is unbroken.

7. **`change_molecule()`** — Add to `ReadOFF`. Implement reference FF loading and
   in-place parameter replacement. This is the primary motivating feature of the merge.

8. **`reference_ff/BLYPSP-4F.off`** — Store the reference water FF in partial `.off`
   format so `change_molecule()` can load it.

9. **OpenMM XML generation** — Port `offtoopenmm/off_to_openmm.py` logic into
   `OpenMMBackend.gen_xml()`, replacing the standalone `OffFileParser` with reads from
   `self._parent.nonbonded`, `self._parent.bonded`, `self._parent.charges`.
   **Note:** `off_to_openmm.py` uses `#define` extraction while `afmtogmx` uses
   `functions.py`-style parsing — verify data-structure equivalence before porting,
   and fill any gaps in `functions.py` rather than working around them.

---

## Full Example Workflow (After Implementation)

```python
import afmtogmx as afm

off = afm.ReadOFF('butanol_water.off')

# Replace refitted water params with the stored standard FF
off.change_molecule(
    mol_name='H2OQM',
    reference_ff='BLYPSP-4F',
    atom_name_map={'OW': 'OW_sp4f', 'HW': 'HW_sp4f'}
)

off.load_charges_from_file('charges.txt')

# GROMACS outputs
off.gmx.set_config(tabpot_prefix='table', scale_C6=True)
off.gmx.gen_nonbonded_tabpot()
off.gmx.write_nonbonded_tabpot()
off.gmx.gen_nonbonded_topology(template_file='template.top', write_to='nonbonded.top')
off.gmx.gen_bonded_topology(template_file='nonbonded.top', write_to='topol.top')

# OpenMM outputs (same parsed data, different output format)
off.openmm.set_config(molname_translations={'H2OQM': 'SOL'})
off.openmm.gen_xml(output='forcefield.xml')
```
