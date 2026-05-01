# Plan: `ReadOFF.populate_bonded` — Bonded Reconstruction From User-Supplied Topology

## Context

For larger molecules built by extending or polymerizing a fitted monomer (e.g., cucurbituril from a glycoluril monomer), the user already knows three things that the package does not:

1. **The new connectivity** — which atom numbers in the larger system are bonded.
2. **The atom typing** — which atom type each new atom carries (reusing the monomer's types so existing nonbonded params transfer).
3. **Which dihedral types were actually fitted** — combinatorial dihedral generation overproduces; only the user's `valid_dihedrals.dat` knows what's real.

A precedent solving the same problem for GROMACS lives in `test/populated_bonded/populate_bonded/`: small `gen_*types.py` scripts walk a bond graph to derive bonds/angles/dihedrals grouped by atom-type signature, then `populate_bonded.py` substitutes parameter values from a `parameters.dat` file. The output is the bonded portion of a GROMACS `.top`. Charges and nonbonded parameters were copied across by hand because GROMACS reuses `[nonbonded_params]` directly.

This plan ports that workflow into afmtogmx as a `ReadOFF` method that returns a **new `ReadOFF` object** representing the larger system. Because the result is a force-field-level transformation (it mutates `bonded`, `nonbonded`, `charges`, `residues`), the operation lives on `ReadOFF` itself — *not* under `off.gmx.` or `off.openmm.`. After the call, the user proceeds with any backend's normal workflow (`new_off.gmx.gen_*`, `new_off.openmm.gen_xml`, etc.).

The atom-removal problem (e.g., dropping the monomer's C3 cap atoms when polymerizing) is handled implicitly: any atom type that disappears from the union of all molecules' atom types in the new object is also dropped from `nonbonded` and `charges`.

---

## Deliverable (user-facing)

```python
import afmtogmx as afm

off = afm.ReadOFF(off_loc='cucurbituril/intra.off')   # monomer + co-fitted water

new_off = off.populate_bonded(
    directory='cb7_topology/',           # contains atoms.dat, bonds.dat,
                                         #          valid_dihedrals.dat, parameters.dat
    new_mol_name='CB7',                  # name for the assembled molecule
    remove_molecules=['C7F'],            # original monomer entries to drop
)

# new_off is a fresh ReadOFF carrying:
#   - 'CB7'    bonded section rebuilt from atoms.dat + bonds.dat + parameters.dat
#   - 'H2OQM'  (and any other co-fitted molecule) copied verbatim from the original
#   - nonbonded filtered to (type, type) pairs whose types still appear somewhere
#   - charges initialized to 0.0 for every type that survives

new_off.load_charges_from_file('cb7_charges.txt')   # standard existing workflow
new_off.gmx.gen_nonbonded_tabpot()
new_off.gmx.gen_bonded_tabpot()
new_off.openmm.gen_xml()
```

The method does one job: **rebuild bonded for the new molecule, keep co-fitted molecules, prune nonbonded to surviving types.** Charges are assigned afterward by the existing `load_charges_from_file` path — no special handling here.

---

## Design

### 1. Input file formats

A single user-supplied directory holds all four files. Atoms and bonds are stored once each (eliminating the three-copy duplication in `test/populated_bonded/populate_bonded/Bonds/`, `/Angles/`, `/Dihedrals/`).

**`atoms.dat`** — atom number → atom type. Atom type must already exist in the original `off`'s atom-type set.
```
  1   C1
  2   C1
  3   N1
  ...
```

**`bonds.dat`** — bonded atom pairs (atom numbers from `atoms.dat`).
```
  1   2
  1   5
  2   3
  ...
```

**`valid_dihedrals.dat`** — dihedral type signatures that should be emitted. One per line. Generated dihedrals not matching any signature (forward or reversed) are silently dropped.
```
C1-C1-N1-C2
N1-C2-N1-C1
O1-C2-N1-C1
```

**`parameters.dat`** — `#define` lines, atom-type signature → interaction subtype + parameters. Underscore count in the signature distinguishes bond/angle/dihedral.
```
#define   C1_C1         HAR   0.15468521   149286.94
#define   C1_C1_H1      HAR   110.04697    384.75972
#define   C1_C1_N1_C2   NCO   0.0          -8.137997   3
```

`HAR` and `NCO` map directly to the existing `bonded[mol]['BON'/'ANG'/'DIH']` interaction-subtype keys; no GROMACS funct-type integers are needed.

### 2. New module: `src/afmtogmx/core/populate_bonded.py`

Module-level helpers, mirroring the `pdb_processing.py` / `xml_generation.py` orchestrator-plus-helpers pattern. No classes — the orchestrator is a single function that the `ReadOFF` method calls.

**Section structure (single file):**
```python
# ---- File parsers ----
def _parse_atoms(path)            -> dict[int, str]              # num -> type
def _parse_bonds(path)            -> list[tuple[int, int]]
def _parse_valid_dihedrals(path)  -> set[str]                    # canonical signatures
def _parse_parameters(path)       -> dict[str, dict]             # by underscore count

# ---- Graph traversal ----
def _build_adjacency(bonds)       -> dict[int, list[int]]
def _find_angles(adjacency)       -> list[tuple[int, int, int]]
def _find_dihedrals(angles, adj)  -> list[tuple[int, int, int, int]]

# ---- Atom-type grouping ----
def _group_bonds_by_type(bonds, atoms_by_num)         -> dict[tuple, list]
def _group_angles_by_type(angles, atoms_by_num)       -> dict[tuple, list]
def _group_dihedrals_by_type(dihedrals, atoms_by_num,
                             valid_signatures)         -> dict[tuple, list]

# ---- Parameter binding ----
def _build_bonded_section(grouped, parameters, kind)  -> dict   # kind = 'BON'/'ANG'/'DIH'

# ---- Exclusions from new bond graph ----
def _generate_exclusions(adjacency, n_atoms)          -> list[list[int]]

# ---- Orchestrator ----
def build_new_molecule_bonded(directory, atom_type_universe) -> dict
    """Read all four files in `directory`, return a fully-populated
    bonded[new_mol_name] dict (ATO/BON/ANG/DIH/EXC), validating that every
    atom type in atoms.dat is present in `atom_type_universe`."""
```

The orchestrator is responsible for (a) calling the parsers, (b) running graph traversal, (c) grouping by type, (d) binding parameters, (e) regenerating exclusions, and (f) returning a single `bonded[mol]` dict shaped exactly like what `gen_md._gen_bonded` produces. It does **not** mutate any `ReadOFF` state — the calling method does that.

### 3. New `ReadOFF` method (in `gen_md.py`)

```python
def populate_bonded(self, directory, new_mol_name, remove_molecules=None):
    """Build a new ReadOFF whose bonded section for `new_mol_name` is
    reconstructed from a user-supplied topology directory.

    Parameters
    ----------
    directory : str
        Path to a directory containing atoms.dat, bonds.dat,
        valid_dihedrals.dat, parameters.dat.
    new_mol_name : str
        Name under which the assembled molecule is stored in the new
        ReadOFF's `bonded`, `charges`, and `residues` dicts.
    remove_molecules : list[str], optional
        Molnames in the current ReadOFF to drop from the new object
        (e.g., the original monomer being replaced). Co-fitted molecules
        not listed here are carried through unchanged.

    Returns
    -------
    ReadOFF
        A new instance. The original `self` is not modified.
    """
```

Implementation sketch (the method body is deliberately short — all heavy lifting is in `populate_bonded.py`):

1. `new = copy.deepcopy(self)` — start from a clone so co-fitted molecules and original `nonbonded` carry across.
2. For each name in `remove_molecules`: drop from `new.bonded`, `new.charges`, `new.residues['Definitions']`, `new.residues['Residues']`.
3. Compute `atom_type_universe = union(types in every surviving molecule's ATO) ∪ types in atoms.dat`.
4. `new.bonded[new_mol_name] = populate_bonded.build_new_molecule_bonded(directory, atom_type_universe)`.
5. Initialize `new.charges[new_mol_name] = {atype: 0.0 for atype in unique_types_in(atoms.dat)}` (mirrors the existing `__init__` charge-initialization convention, which keys by atom type via `pair[1]`).
6. Recompute `surviving_types = union(types across all molecules in new.bonded)`.
7. `new.nonbonded = {pair: params for pair, params in new.nonbonded.items() if pair[0] in surviving_types and pair[1] in surviving_types}`.
8. Add `new_mol_name` to `new.residues` using the same one-liner pattern as `__init__`.
9. **Rewire backends** so they reference the new instance: `new.gmx = GROMACSBackend(new); new.openmm = OpenMMBackend(new)` (deepcopy preserves backend objects but they'd point at the *clone's* attributes via `self._parent` — verify during implementation; if they hold a direct reference to the original `ReadOFF`, replace them outright).
10. `return new`.

### 4. Bonded-section construction details

For each grouped tuple `(type_signature, [atom_tuple, atom_tuple, ...])`:
- Look up the parameter entry whose name matches the signature (try forward and reversed, mirroring existing `populate_bonded.py` logic at lines 100–106 / 163–168 / 226–230).
- The param entry's subtype string (`HAR`/`NCO`) becomes the inner key under `BON`/`ANG`/`DIH`.
- The param values become the parameter tuple used as the dict key. The list of atom-number tuples becomes the value.

Result shape matches the existing data model exactly:
```python
new.bonded['CB7']['BON']['HAR'] = {(0.15469, 149286.94): [[1, 2], [12, 13], ...], ...}
new.bonded['CB7']['ANG']['HAR'] = {(110.05, 384.76):     [[1, 2, 11], ...], ...}
new.bonded['CB7']['DIH']['NCO'] = {(0.0, -8.138, 3):     [[1, 2, 3, 84], ...], ...}
```

**Fail-fast on any missing parameter — no silent fallthroughs anywhere.** Producing a force field with incomplete bonded coverage would silently bias a simulation; the cost of a verbose error is far below the cost of a wrong production run. Concretely:

- **Every generated bond signature** must resolve to a `parameters.dat` entry (forward or reversed). If not → collect into a missing-set; do not emit a partial bonded section.
- **Every generated angle signature** must resolve. Same treatment.
- **Every dihedral signature listed in `valid_dihedrals.dat`** must resolve. (Dihedrals generated by graph traversal but absent from `valid_dihedrals.dat` are filtered out — that's the filter's purpose, not an error. But if the user *did* list a signature in `valid_dihedrals.dat` and forgot the matching `#define`, that's an error.)
- **Every atom type in `atoms.dat`** must already exist in `atom_type_universe` (the union of types across the original `off`'s surviving molecules). A new type would have no nonbonded params, so we refuse rather than silently emit a force field with missing pair interactions.

After processing every category, if the missing-set is non-empty, raise `ValueError` once with the **full categorized list** of missing items (missing bond signatures, missing angle signatures, missing dihedral signatures, unknown atom types). One error message, all gaps named — the user fixes them in one editing pass instead of running the method N times to discover N missing entries.

### 5. ATO format and atom naming

`atoms.dat` has only `(number, type)` — no atom name column. Following the existing `__init__` convention where charges are keyed by atom *type* (`pair[1]`), we use `atom_name = atom_type` for every entry:

```python
new.bonded['CB7']['ATO']['All'] = {1: ('C1', 'C1'), 2: ('C1', 'C1'), 3: ('N1', 'N1'), ...}
new.bonded['CB7']['ATO']['Virtual'] = {}     # no virtuals from this path
```

Multiple atoms sharing a name is fine because ATO is keyed by atom number; the existing charges dict already collapses by type.

### 6. NETF / TORQ handling

The new molecule has none — `atoms.dat` is the user's complete atom list. Skip the NETF/TORQ append step entirely for `new_mol_name`. Co-fitted molecules retain whatever NETF/TORQ entries they had in the original.

---

## Files to Create / Modify

| File | Action | Purpose |
|------|--------|---------|
| `src/afmtogmx/core/populate_bonded.py` | **CREATE** | Parsers, graph traversal, type grouping, parameter binding, exclusion generation |
| `src/afmtogmx/core/gen_md.py` | **MODIFY** | Add `ReadOFF.populate_bonded()` method (orchestrates the above) |
| `test/populated_bonded/populate_bonded/...` | **REFERENCE ONLY** | Existing scripts stay as-is — they are the reference implementation we are porting |
| `CLAUDE.md` | **MODIFY** | Document the new method + input file formats |

### Existing code reused (no duplication)

- **`copy.deepcopy`** — clone the original `ReadOFF` so co-fitted molecules and the existing nonbonded dict come along for free.
- **`functions._remove_netf_torq_atname` / `_atnum`** (`functions.py:907`/`929`) — used when building the residues entry for the new molecule, matching the existing one-liner in `gen_md.py:76`.
- **Existing graph-traversal logic** in `test/populated_bonded/populate_bonded/Angles/gen_angletypes.py` (`find_angles`) and `Dihedrals/gen_dihedraltypes.py` (`find_dihedrals` + `valid_dihedrals` filter) — ported directly into module-level helpers.
- **Existing parameter-format convention** in `populate_bonded.py:6-43` — `#define <name> <type> <values...>` with underscore count distinguishing bond/angle/dihedral. Identical parser, just returns a dict instead of writing a `.top`.
- **Existing canonical-signature normalization** in `populate_bonded.py:45-63` — alphabetical sort for bonds, endpoint-swap canonicalization for angles, reverse-canonicalization for dihedrals.

### Architectural placement

Method lives on `ReadOFF`, exactly like `change_molecule` (`gen_md.py:325`) and `load_charges_from_file`. It is **not** placed on `off.gmx` or `off.openmm` because the operation is force-field-level (it mutates the same dicts that backends consume) and is independent of any output format.

---

## Verification

1. **Parser round-trip** — write a synthetic `atoms.dat`/`bonds.dat`/`valid_dihedrals.dat`/`parameters.dat` for a small toy system; assert each `_parse_*` helper returns the expected structure.
2. **Graph traversal against the cucurbituril reference** — feed `test/populated_bonded/populate_bonded/Bonds/atoms.dat` + `bonds.dat` into `_find_angles` and `_find_dihedrals`, confirm the count of unique angle/dihedral type signatures matches what `gen_angletypes.py` / `gen_dihedraltypes.py` produce in the existing `angletypes.dat` / `dihedraltypes.dat`.
3. **Parameter binding** — run `_build_bonded_section` against `parameters.dat`, confirm the resulting dict matches the parameter values printed in `bonded.top`.
4. **Type-removal cleanup** — supply an `atoms.dat` lacking type `C3`, confirm every `('C3', *)` and `(*, 'C3')` key is absent from `new.nonbonded` and that `'C3'` is absent from `new.charges['CB7']`.
5. **Co-fitted pass-through** — confirm a co-fitted molecule (e.g., `H2OQM`) and its cross-pair nonbonded entries (e.g., `('OW', 'C1')`) survive untouched in `new_off`.
6. **Missing-parameter fail-fast** — delete one `#define` line from `parameters.dat`, confirm `populate_bonded()` raises with the exact missing signature in the message.
7. **Backend round-trip** — call `new_off.gmx.gen_nonbonded_tabpot()` and `new_off.gmx.gen_bonded_tabpot()` on a populated CB7 example; confirm the generated tables include only surviving types and that the bonded `.top` matches the reference `bonded.top` line-for-line in counts (parameter values follow `parameters.dat`).
8. **Existing regression tests** — `pytest test/` continues to pass; this plan adds files but does not change existing behavior.

---

## Explicitly Out of Scope (defer to follow-up plans)

- **Charges from a by-type `charges.dat`**. The new `ReadOFF` initializes charges to 0.0 per type; the user calls the existing `load_charges_from_file` afterward. A by-type loader can be added later if needed but is not part of this method.
- **Multiple new molecules in one call**. One `populate_bonded` invocation = one `(directory, new_mol_name)`. Building two new molecules means two calls (chainable since each returns a new `ReadOFF`).
- **Heuristic dihedral inference** without `valid_dihedrals.dat`. The file is required; we do not try to guess which combinatorial dihedrals were actually fitted.
- **Polymer/ring assembly conveniences** (the `repeat()` API in `PLANNING_MOLECULE_BUILDER.md`). The `bonds.dat` here is the user's responsibility to construct; this plan is the parameter-binding layer only.
- **SMILES-driven discovery** (`PLANNING_SMILES_BANK.md`). This plan assumes the user has already chosen which monomer's parameters to extend and how to connect units; SMILES match would sit upstream of `populate_bonded`.
- **Bonded-term fitting**. We bind precomputed parameters; we do not generate or refit them.

---

## Relationship to the other planning docs

The three plans address different layers of the same broad problem (build a new force field from existing fitted ones):

- **`PLANNING_POPULATE_BONDED.md`** (this plan) = **lowest layer, highest priority.** Mechanical port of a working precedent. Requires the user to supply atom typing, bond list, dihedral filter, and parameter values explicitly. Smallest scope, fastest to ship, immediately useful for the cucurbituril case.
- **`PLANNING_MOLECULE_BUILDER.md`** = **assembly layer above this one.** Atom-number-level fragment selection, junction-term auto-generation, polymer `repeat()` convenience. Could be built on top of `populate_bonded` once that exists, with `repeat()` emitting the `atoms.dat`/`bonds.dat` that `populate_bonded` consumes.
- **`PLANNING_SMILES_BANK.md`** = **discovery layer above both.** SMILES-driven substructure search across a bank of fitted systems, system-context coverage tiers, gap reporting. The most ambitious but also the most front-loaded. Could feed selections into either of the layers below.

Recommended ordering: ship `populate_bonded` first (small, validates against an existing GROMACS-side workflow), see whether it covers enough use cases on its own before committing to the larger assembly/discovery layers.
