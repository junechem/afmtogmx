# Plan: MoleculeBuilder — Fragment Assembly System

## Context

The user needs to build new molecules from fragments of existing fitted force fields. Two concrete use cases:

1. **Polymer assembly**: Take a fitted monomer (e.g., C7F glycoluril from cucurbituril.off), repeat it N times to form a polymer or ring. Remove capping atoms, add inter-unit bonds, borrow missing junction parameters.

2. **Molecule extension**: Take a fitted molecule (e.g., 1-butanol), extend it by grafting fragments from other sources (e.g., CH2 groups) with some parameters borrowed from a different force field.

Both are the same problem: **select atoms from existing molecules, combine them, add new bonds, supply missing parameters, produce a new force-field-ready object.**

Charges are the user's responsibility — the tool warns about neutrality but doesn't solve it.

---

## Architecture: Three Components

### 1. `MoleculeData` — Result container

**File**: `src/afmtogmx/core/molecule_data.py` (new)

Duck-type compatible with `ReadOFF`. Holds the same four data dicts and wires up both backends:

```python
class MoleculeData:
    def __init__(self, mol_name, bonded, nonbonded, charges):
        self.bonded = bonded          # {mol_name: bonded_dict}
        self.nonbonded = nonbonded    # {(type1, type2): {int_type: [params]}}
        self.charges = charges        # {mol_name: {atname: charge}}
        self.residues = {             # auto-built from bonded
            "Definitions": {...},
            "Residues": {...}
        }
        self.gmx = GROMACSBackend(self)
        self.openmm = OpenMMBackend(self)

    def load_charges_from_file(self, file_path):
        # Same logic as ReadOFF.load_charges_from_file
```

Backends access `self._parent.bonded` etc. — `MoleculeData` works as drop-in parent. **No backend changes needed.**

### 2. `Fragment` — Internal atom subset

**File**: `src/afmtogmx/core/molecule_builder.py` (new)

Created from source `ReadOFF` + molecule name + atoms to keep/remove. Responsibilities:

- **Atom selection** by atom number (unambiguous — names repeat)
- **Bonded term filtering**: keep only terms where ALL referenced atoms are in the subset
- **Renumbering**: old → new contiguous integers, stored in `renumber_map`
- **NETF/TORQ stripping** (re-added at final build)
- **Virtual site handling**: drop if any referenced atom removed

```python
class Fragment:
    def __init__(self, source, mol_name, keep_atoms=None, remove_atoms=None, type_map=None):
        # Exactly one of keep_atoms/remove_atoms must be provided (lists of int)
        # type_map: optional {old_type: new_type} for cross-source assembly

    # After construction, Fragment holds:
    self.atoms         # {int: (atname, attype)} — renumbered
    self.virtuals      # virtual site defs, renumbered
    self.bonds         # {(params): [[at1, at2], ...]} — renumbered, per subtype
    self.angles        # same pattern
    self.dihedrals     # same pattern
    self.bd3           # same pattern
    self.cdihedrals    # same pattern
    self.renumber_map  # {old_int: new_int}
    self.reverse_map   # {new_int: old_int}
    self.source_mol    # original molecule name
    self.nonbonded     # reference to source nonbonded (for type-pair lookup later)
```

**Bonded term filtering algorithm**: For each section (BON/ANG/DIH/BD3/CDI), for each subtype (HAR/QUA/NCO/etc.), for each parameter tuple -> filter atom lists to keep only entries where every atom is in the kept set. Remove empty parameter tuples. Remove empty subtypes.

**ATO key normalization**: ATO uses string keys in ReadOFF (`'1'`, `'2'`), but BON/ANG/DIH use int atom numbers. Fragment normalizes everything to int internally. ATO keys converted back to strings in final output.

### 3. `MoleculeBuilder` — User-facing builder

**File**: `src/afmtogmx/core/molecule_builder.py` (same file)

---

## Full API Design

```python
class MoleculeBuilder:
    def __init__(self, mol_name):
        """Create a builder for a new molecule.
        mol_name: name for the assembled molecule in output dicts
        """

    # -- Fragment Management --

    def add_fragment(self, source, mol_name, keep_atoms=None, remove_atoms=None,
                     label=None, type_map=None):
        """Add a fragment from a ReadOFF or MoleculeData source.

        source: ReadOFF or MoleculeData object
        mol_name: molecule name within the source
        keep_atoms / remove_atoms: list of int atom numbers (exactly one required)
        label: string label for this fragment (auto-generated if None: 'frag_0', 'frag_1', ...)
        type_map: optional {old_type: new_type} dict for renaming atom types

        Returns: the label (str)
        """

    # -- Bond Declaration --

    def add_bond(self, frag1_label, atom1, frag2_label, atom2):
        """Declare a new bond between fragments.

        atom1/atom2: atom numbers in the ORIGINAL source molecule (before renumbering).
        The builder resolves these through the fragment's renumber_map.
        """

    # -- Parameter Sources --

    def set_params(self, source, mol_name=None):
        """Register a parameter source for junction term lookup.

        source: ReadOFF or MoleculeData -- its bonded dict is searched for
                matching atom-type patterns when junction terms need parameters.
        mol_name: if specified, only search this molecule's parameters.
                  If None, search all molecules in the source.

        Can be called multiple times. Sources are searched in registration order.
        """

    def add_bond_params(self, bond_type, atom_types, params):
        """Explicitly add bond parameters.
        bond_type: 'HAR' or 'QUA'
        atom_types: tuple of atom type names, e.g. ('N1', 'C3')
        params: tuple of parameter values, e.g. (1.47, 600.0)
        """

    def add_angle_params(self, angle_type, atom_types, params):
        """Explicitly add angle parameters.
        angle_type: 'HAR' or 'QUA'
        atom_types: tuple of 3 atom type names
        params: parameter tuple
        """

    def add_dihedral_params(self, dih_type, atom_types, params):
        """Explicitly add dihedral parameters.
        dih_type: 'HAR', 'NCO', or 'COS'
        atom_types: tuple of 4 atom type names
        params: parameter tuple
        """

    # -- Inspection --

    def show_atoms(self, source=None, mol_name=None):
        """Print atom table. Two modes:

        1. show_atoms(source=off, mol_name='C7F')
           -- Show atoms from a source molecule (for planning which atoms to keep/remove)

        2. show_atoms()
           -- Show all atoms across all added fragments, with columns:
             Fragment | Orig# | New# | Name | Type
        """

    def show_bonds(self):
        """Print existing bonds (from fragments) and declared new bonds."""

    def show_missing_params(self):
        """Dry-run parameter lookup. Print which junction terms have params
        and which are still missing. Helps the user know what to supply."""

    # -- Polymer Convenience --

    def repeat(self, source, mol_name, n_copies,
               remove_atoms_per_copy=None,
               remove_atoms_first=None,
               remove_atoms_last=None,
               link_atoms=None,
               topology='chain'):
        """Convenience for polymer/ring assembly.

        Creates n_copies fragments, each from source/mol_name.

        remove_atoms_per_copy: atoms removed from ALL copies (the standard caps)
        remove_atoms_first: ADDITIONAL atoms removed only from the first copy
        remove_atoms_last: ADDITIONAL atoms removed only from the last copy
        link_atoms: list of (atom_in_prev, atom_in_next) tuples defining
                    how adjacent copies connect. Atom numbers are original.
        topology: 'chain' or 'ring'
            - 'chain': first copy uses remove_atoms_first, last uses remove_atoms_last
            - 'ring': all copies use remove_atoms_per_copy only; last links back to first

        Fragment labels are auto-generated as 'unit_0', 'unit_1', etc.

        Note: link_atoms connect atoms that are KEPT (not removed). The atoms
        that connect may be the same number in each copy (e.g., N atoms at
        positions 3 and 5) -- the fragment label disambiguates them.
        """

    # -- Build --

    def build(self, auto_junctions=True):
        """Assemble everything into a MoleculeData object.

        auto_junctions: if True, automatically generate angle and dihedral
                        terms at junction bonds by examining the bond graph
                        neighbors on each side. If False, only explicitly added
                        terms are used at junctions.

        Steps:
        1. Concatenate fragment atom tables (renumbering sequentially)
        2. Concatenate all bonded terms with renumbered atom references
        3. Add declared inter-fragment bonds
        4. If auto_junctions: generate angle/dihedral terms at junctions
        5. Look up parameters for all junction terms
        6. Regenerate exclusion list from bond graph (1-2 and 1-3 neighbors)
        7. Append NETF + TORQ atoms
        8. Collect nonbonded terms (union of all atom-type pairs from sources)
        9. Build charges dict (default 0.0, warn about neutrality)
        10. Build residues dict
        11. Return MoleculeData

        Raises ValueError if any junction terms are missing parameters.
        The error message lists exactly which terms need params and their types.
        """
```

---

## Parameter Lookup Algorithm

When `build()` encounters a junction bond between atom types `(T1, T2)`:

1. **Check explicit params** -- `add_bond_params('HAR', ('T1', 'T2'), ...)` matches directly.

2. **Search registered sources** -- For each source (in order):
   - Scan all molecules' `BON['HAR']` entries
   - For each parameter tuple, check if any atom pair `[a, b]` has types matching `(T1, T2)` or `(T2, T1)`
   - If found, use those parameters
   - First match wins

3. **If not found** -- Record as missing. After all junction terms are processed, raise one error listing all missing terms.

Same algorithm for angles (match 3-type patterns) and dihedrals (match 4-type patterns).

**Type matching for angles/dihedrals**: The atom-type sequence must match either forward or reversed. For angles `(T1, T2, T3)`, match against `(T1, T2, T3)` or `(T3, T2, T1)`. For dihedrals, match forward or reversed.

---

## Junction Term Auto-Generation

When `auto_junctions=True` and a new bond A-B is declared:

1. **Build bond graph** from all fragment bonds + declared bonds
2. **For angles**: For each neighbor A' of A (where A'!=B), create angle term A'-A-B. For each neighbor B' of B (where B'!=A), create angle term A-B-B'.
3. **For dihedrals**: For each A'-A-B, for each B' of B (B'!=A), create dihedral A'-A-B-B'. For each A-B-B', for each A' of A (A'!=B), create dihedral A'-A-B-B'. Also extend one more step for proper 4-atom chains.

Each generated term is looked up using the parameter lookup algorithm above.

---

## Exclusion Regeneration Algorithm

```python
def _generate_exclusions(all_bonds, n_atoms):
    """Build EXC list from bond graph.

    all_bonds: list of [at1, at2] pairs (ints)
    n_atoms: total real atoms (excluding NETF/TORQ)

    Returns: list of lists in .off EXC format (upper-triangle)
    """
    # Build adjacency
    neighbors = defaultdict(set)
    for a1, a2 in all_bonds:
        neighbors[a1].add(a2)
        neighbors[a2].add(a1)

    # For each atom, collect 1-2 and 1-3 neighbors with higher numbers
    exclusions = []
    for atom in range(1, n_atoms + 1):
        excluded = set()
        for n1 in neighbors[atom]:
            if n1 > atom:
                excluded.add(n1)
            for n2 in neighbors[n1]:
                if n2 > atom and n2 != atom:
                    excluded.add(n2)
        if excluded:
            exclusions.append([atom] + sorted(excluded))
    return exclusions
```

**Note**: Also add NETF/TORQ to every atom's exclusion list (matching existing .off convention).

---

## Nonbonded Collection

At build time:
1. Collect all unique atom types across all fragments
2. For each unique pair `(T1, T2)` (sorted):
   - Search all fragment sources for nonbonded params for this pair
   - If found in multiple sources with DIFFERENT params -> error, tell user to use `type_map`
   - If found -> include in output
   - If not found -> warn (might be intentional if types don't interact)

---

## Charge Handling

- Initialize all charges to 0.0
- User sets charges via `result.charges[mol_name][atname] = value` or `result.load_charges_from_file()`
- At `build()` time: print total charge. If not zero (within tolerance), print warning:
  `"Warning: Total charge is {x:.4f}. Adjust charges for neutrality."`

---

## Files to Create/Modify

| File | Action | Purpose |
|------|--------|---------|
| `src/afmtogmx/core/molecule_data.py` | CREATE | MoleculeData container |
| `src/afmtogmx/core/molecule_builder.py` | CREATE | Fragment + MoleculeBuilder |
| `src/afmtogmx/__init__.py` | MODIFY | Export MoleculeBuilder, MoleculeData |

**No changes to existing modules**: gen_md.py, gmx_backend.py, openmm_backend.py, functions.py, topology.py, tabulated_potentials.py, residues.py

### Existing functions to reuse
- `functions.gen_empty_bonded()` (`functions.py:175`) -- bonded dict template
- `functions.gen_empty_nonbonded()` (`functions.py:237`) -- nonbonded dict template
- `functions._remove_netf_torq_atname()` (`functions.py:907`) -- filter special atoms by name
- `functions._remove_netf_torq_atnum()` (`functions.py:929`) -- filter special atoms by number
- `functions._filter_bonded()` (`functions.py:844`) -- remove empty sections

---

## Implementation Phases

### Phase 1: Core data flow
1. Create `MoleculeData` with backend wiring + `load_charges_from_file()`
2. Create `Fragment` -- atom selection, bonded term filtering, renumbering
3. Create `MoleculeBuilder` with `add_fragment()`, `add_bond()`, basic `build()`
   - Build concatenates fragments, adds declared bonds, generates EXC, collects nonbonded
   - No parameter lookup yet -- junction bonds must have params from fragment sources
4. Test: parse ethane -> split into two CH3 fragments -> reassemble -> verify bonded dict

### Phase 2: Parameter management
5. Add `set_params()` and explicit param methods (`add_bond_params`, etc.)
6. Add parameter lookup algorithm in `build()`
7. Add junction term auto-generation (`auto_junctions=True`)
8. Add `show_missing_params()` for dry-run diagnostics
9. Test: molecule extension scenario

### Phase 3: Polymer convenience + display
10. Add `repeat()` method with chain/ring topology
11. Add `show_atoms()`, `show_bonds()` display methods
12. Add `type_map` support on `add_fragment()`
13. Add charge neutrality warning
14. Test: build 2-unit and 7-unit cucurbituril

### Phase 4: Final integration
15. Export from `__init__.py`
16. Write regression tests
17. Update CLAUDE.md with MoleculeBuilder docs

---

## Verification

1. **Roundtrip test**: Parse ethane.off -> split into two CH3 -> reassemble -> compare bonded/nonbonded with original
2. **Fragment test**: Create fragment from C7F with removed atoms, verify bonded terms filtered correctly
3. **Renumbering test**: Verify all atom references in BON/ANG/DIH update correctly
4. **EXC test**: Verify regenerated exclusions match expected pattern
5. **Parameter lookup test**: Register a source, verify junction terms find correct params
6. **Integration test**: Build 2-unit chain from C7F -> verify `result.gmx.gen_nonbonded_tabpot()` works
7. **Ring test**: Build 7-unit cucurbituril ring -> verify ring closure bonds present

---

## Edge Cases to Handle

- **String vs int ATO keys**: Normalize to int internally, convert back on output
- **Virtual sites**: Carry through if all referenced atoms kept; drop otherwise
- **BD3/CDI sections**: Filter with same logic as BON/ANG but with 3/5 atom refs
- **Empty parameter subtypes**: Clean up after filtering (reuse `_filter_bonded()`)
- **Duplicate atom types across sources**: `type_map` renames to avoid collisions
- **NETF/TORQ**: Strip from fragments, add single pair to final molecule

---

## Internal Data Structure Reference

```python
# ATO: string keys -> (atname, attype) tuples
bonded['MolName']['ATO']['All'] = {'1': ('C1', 'C1'), '2': ('H1', 'H1')}
bonded['MolName']['ATO']['Virtual'] = {(5, 'O', 'O1'): ('def1', 'def2')}

# BON: param tuples -> lists of atom-number pairs (ints)
bonded['MolName']['BON']['HAR'] = {(0.123, 456.7): [[1, 2], [3, 4]]}
bonded['MolName']['BON']['QUA'] = {(0.1, 0.2, 0.3, 0.4): [[5, 6]]}

# ANG: param tuples -> lists of atom-number triples (ints)
bonded['MolName']['ANG']['HAR'] = {(110.0, 50.0): [[1, 2, 3]]}

# DIH: param tuples -> lists of atom-number quads (ints)
bonded['MolName']['DIH']['NCO'] = {(10.0, 0.5, 0.6): [[1, 2, 3, 4]]}

# EXC: list of lists (upper-triangle)
bonded['MolName']['EXC'] = [[1, 2, 3, 4], [2, 3, 4], [3, 4]]

# Nonbonded: sorted type-pair tuples -> interaction dicts
nonbonded[('C1', 'H1')] = {'EXP': [[params]], 'SRD': [[params]], ...}

# Charges: molname -> atname -> float
charges['MolName'] = {'C1': 0.0, 'H1': 0.0}

# Residues
residues = {
    "Definitions": {"MolName": {"All": {1: ('C1', 'C1'), ...}}},  # netf/torq filtered
    "Residues": {"MolName": {"All": [[1, 2, 3, ...]]}}            # atom numbers
}
```
