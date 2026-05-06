# Plan: SMILES-Aware Force Field Bank — Foundational Layer

## Context

The user wants to build new molecular force fields by reusing parameters from a bank of previously fitted ones, driven by SMILES rather than manual atom-number selection (the approach in the companion `PLANNING_MOLECULE_BUILDER.md`). Two physics-driven constraints shape the design:

1. **Each .off file IS a fitted system, not a single molecule.** Inspection of `test/sample_off_files/butanediol_intra.off` confirms its `nonbonded` dict contains cross-molecule pairs like `('OW', 'C1')` and `('HW', 'O1')` — water-butanol terms fitted *together*. So the bank's natural unit is "fitted system," and pairwise transferability depends on whether the *system context* matches.

2. **Charge neutrality is unenforced today.** `ReadOFF.charges` defaults all atoms to 0.0 with no validation; any builder must surface charge totals so the user can fix imbalances before simulation.

This plan delivers the **foundational discovery layer**: SMILES annotation format, bank indexing, substructure search across systems, and a coverage-report tool that tells the user *given my target system, what does my bank cover and what are the gaps?* — without yet generating output force-field files. Once the user trusts the search results against real .off files, a follow-up plan layers the actual FF assembly on top, reusing the mechanics sketched in `PLANNING_MOLECULE_BUILDER.md`.

---

## Deliverable (user-facing)

```python
import afmtogmx as afm

bank = afm.load_bank('path/to/off_files/')

report = bank.coverage_report(
    target_system=[
        ('CCCC(O)O', 'BUOH', 1),     # 1 butanediol
        ('O',         'WAT', 343),    # 343 waters
    ],
)
print(report)
# Per molecule:  best source .off, substructure alignment, missing atoms
# Per pair:      best source .off, system-context tier (1/2/none), source atom-type pair
# Per molecule:  total charge (× count), neutrality flag
# System total:  net charge, neutrality flag
# Gaps:          explicit list of pairs/atoms with no acceptable source
```

The user iterates: annotate more .off files → rerun the report → resolve gaps. **No force-field files are written by this plan.**

---

## Design

### 1. SMILES annotations live in .off header comments

Per-molecule SMILES + per-atom mapping to RDKit canonical indices. Hydrogens are explicit (`Chem.AddHs`).

```
# SMILES: BUUQM CCCC(O)O
# ATOMMAP: BUUQM C1=0 C2=1 C3=2 C4=3 O1=4 H1=5 H2=6 ...
# SMILES: H2OQM O
# ATOMMAP: H2OQM O=0 H1=1 H2=2

[FF_INPUT]
...
```

`#`-prefixed lines are already ignored by the existing parser (verified during exploration), so this is purely additive — no changes to `gen_md.py` or `functions.py` required.

### 2. New module: `src/afmtogmx/core/smiles_bank.py`

Stateless module-level helpers + lightweight dataclasses, mirroring the orchestrator pattern used in `xml_generation.py` / `pdb_processing.py`.

**Data model**:
```python
@dataclass
class MolEntry:
    name: str                    # .off molecule name (e.g., 'BUUQM')
    smiles: str
    rdmol: Chem.Mol              # with explicit Hs
    atom_map: dict[str, int]     # off_atom_name → rdkit atom idx
    atom_types: dict[str, str]   # off_atom_name → off atom_type
    charges: dict[str, float]    # off_atom_name → charge (0.0 default)

@dataclass
class BankEntry:
    path: str                          # path to .off file
    readoff: ReadOFF                   # the parsed force field
    molecules: dict[str, MolEntry]     # molname → MolEntry
    nonbonded_pairs: set[tuple[str,str]]  # all (atype, atype) pairs in this system

@dataclass
class CoverageReport:
    molecule_matches:  list[MoleculeMatch]   # one per target molecule
    pair_matches:      list[PairMatch]       # one per needed nonbonded pair
    charge_summary:    ChargeSummary
    gaps:              list[Gap]             # everything that fails to match

    def __str__(self) -> str: ...   # formatted text report
```

**Public API**:
- `parse_off_smiles_metadata(off_path)` → `dict[molname, (smiles, atom_map)]`. Reads the `# SMILES:` and `# ATOMMAP:` header comments before any `[SECTION]` line. Raises if a SMILES line lacks a matching ATOMMAP for any molecule.
- `load_bank(directory)` → `Bank`. Walks the directory for `*.off`, parses metadata, instantiates `ReadOFF` per file, validates that every annotated atom appears in the .off ATO section.
- `Bank.coverage_report(target_system)` → `CoverageReport` (see algorithm below).

**Internal helpers** (each kept small for unit-testability):
- `_substructure_align(target_rdmol, source_rdmol)` → list of `{target_idx: source_idx}` mappings via RDKit's `GetSubstructMatches`
- `_translate_alignment_to_atom_names(alignment, target_atom_map, source_mol_entry)` → maps the rdkit-index alignment back to .off atom names on each side
- `_match_molecule(target_smiles, bank)` → ranked candidates: exact > superstructure > substructure-only
- `_match_pair(target_pair_atoms, mol_alignments, bank)` → best source for one pair, with system-context tier (defined below)
- `_compute_charges(target_system, mol_alignments)` → per-mol totals + system total

### 3. Search algorithm

For each target molecule:
1. Parse SMILES, add explicit Hs.
2. Iterate bank entries, for each `MolEntry`:
   - Try `target.GetSubstructMatch(source.rdmol)` (source ⊆ target — partial coverage)
   - Try `source.rdmol.GetSubstructMatch(target)` (target ⊆ source — full coverage from a larger fitted molecule)
   - Score: full coverage > superstructure-coverage > partial.
3. Return the ranked list. Top match supplies the per-atom alignment used downstream.

For each needed nonbonded pair (intra-A, intra-B, **and** every A-B cross combination):
1. Identify the two atoms by their `(target_molecule, target_atom_idx)`.
2. **Tier 1 — system match**: bank entry where *both* contributing molecules align AND the corresponding source atom-type pair exists in that entry's `nonbonded_pairs`. These params were fitted in matching context — high confidence.
3. **Tier 2 — cross-context match**: each atom found in some bank entry, but not the *same* entry. Use a tier-2 source and flag the report.
4. **Tier 3 — gap**: log it; the report's `gaps` list collects all of these. Per the user's chosen behavior, the report makes the gap loud but does not raise — the report itself is the artifact.

### 4. Charge handling

Per-target-molecule charge = sum of (mapped source atom charges × molecule count). System total = sum across all target molecules. Report emits:
- Per-molecule total + neutrality flag (`|charge| < 1e-4`)
- System total + neutrality flag
- A warning if any molecule charge is non-integer (suggests fitting issue or missing charge file)

No automatic adjustment — the user owns charge resolution.

### 5. Annotation tooling (minimum viable)

A single helper to make the one-time annotation manageable:
- `bank.suggest_atommap(off_path, mol_name, smiles)` → prints a candidate mapping by element-and-degree heuristics, which the user verifies and pastes into the .off header. Not fully automatic — atom-naming conventions in .off are user-driven and shouldn't be guessed silently.

---

## Files to Create / Modify

| File | Action | Purpose |
|------|--------|---------|
| `src/afmtogmx/core/smiles_bank.py` | **CREATE** | Bank, entries, search, coverage report |
| `src/afmtogmx/__init__.py` | **MODIFY** | Export `load_bank` |
| `test/sample_off_files/water_intra.off` | **MODIFY** | Add SMILES/ATOMMAP header comments (annotate H20QM, water variants) |
| `test/sample_off_files/butanediol_intra.off` | **MODIFY** | Add SMILES/ATOMMAP header comments (BUUQM, H2OMM, H2OQM) |
| `pyproject.toml` | **MODIFY** | Add `rdkit` dependency |
| `CLAUDE.md` | **MODIFY** | Document SMILES annotation format and `load_bank` / `coverage_report` API |

### Existing code reused (no duplication)

- **`ReadOFF.__init__`** (`src/afmtogmx/core/gen_md.py`) — each bank entry holds a parsed `ReadOFF`; we don't re-implement .off parsing.
- **`ReadOFF.nonbonded`** dict structure (sorted atom-type tuple keys, confirmed at `functions.py:835`) — directly indexed by the search algorithm.
- **`ReadOFF.charges`** (`gen_md.py:72-74`) — directly read for the charge summary.
- **`change_molecule`** (`gen_md.py:325-482`) — *not* called by this plan, but its swap-with-rename logic is the proven precedent that the follow-up FF-generation plan will build on.

---

## Verification

1. **Annotate two sample files** — manually add `# SMILES:` and `# ATOMMAP:` lines to `water_intra.off` and `butanediol_intra.off`. Run `parse_off_smiles_metadata` and confirm the maps round-trip.
2. **Bank load** — `bank = load_bank('test/sample_off_files/')`. Verify both files are picked up and every annotated atom name is present in each entry's ATO section.
3. **Exact-system coverage** — target = 1 BUUQM + 343 H2OQM (the system that butanediol_intra.off was fitted for). Report should show 100% coverage, all pairs at tier 1, charges as currently set in the file.
4. **Partial substructure coverage** — target = propanol + water. Propanol must match as a substructure of butanediol; report should show partial alignment with the missing CH2 listed under per-molecule gaps; intra-propanol pairs available, propanol-water cross pairs at tier 1 where covered.
5. **No-match case** — target = methane + water. Methane has no source; report must list per-pair gaps with the closest-source suggestion (e.g., "no source contains a CH4-equivalent fragment").
6. **Charge totals** — populate non-zero charges via `load_charges_from_file` on a bank file, rerun the report, confirm system total reflects per-molecule × count.
7. **Existing regression tests** — `pytest test/` must still pass (this plan is purely additive; no existing module is changed).

---

## Explicitly Out of Scope (defer to next plan)

- Writing assembled `.off` / GROMACS / OpenMM force-field files from a coverage report.
- Fragment renumbering, junction-term auto-generation, exclusion regeneration — the mechanics in `PLANNING_MOLECULE_BUILDER.md` will be revisited when building the FF-emission layer beneath this discovery layer.
- Automatic charge neutralization or fitting suggestions.
- Auto-generating SMILES annotations from existing .off files (heuristic-only `suggest_atommap` is the most we attempt).
- Polymer/ring `repeat()` convenience — the SMILES-level expression of polymers is a separate design problem (user provides repeated-unit SMILES vs. polymer SMARTS — needs its own discussion).

---

## Relationship to `PLANNING_MOLECULE_BUILDER.md`

The two plans are **complementary, not competing**:

- This plan (SMILES bank) = **discovery layer** — chemically intuitive UX, system-aware coverage analysis. The user's first interaction with the tool.
- `PLANNING_MOLECULE_BUILDER.md` = **assembly layer** — atom-renumbering, junction-term generation, exclusion regeneration, residue-dict construction. Required mechanics regardless of whether the user picks atoms manually or via SMILES match.

The intended ordering: ship this foundational layer first → user validates that the search/coverage logic produces sensible matches against their real bank → next plan builds the FF emission, driven by coverage-report results, leaning on the assembly mechanics from `PLANNING_MOLECULE_BUILDER.md`. The manual atom-selection API from that older plan may survive as a low-level fallback for cases the SMILES path can't handle (e.g., novel chemistry with no bank coverage).
