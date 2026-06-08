# `functions.py` — `.off` parsing utilities + interaction math

Low-level helpers used mainly by `ReadOFF.__init__` to parse the `.off` file, plus
the pairwise interaction functions used during tabulated-potential generation.

**Module global:** `total_bonded_added` (int) — running counter incremented while
parsing fitted bonded sections, used to index into the unsorted intra-potential
list. `ReadOFF._gen_bonded()` resets it to 0 after each molecule pass.

## Interaction functions (potential + force)
All return `(potential, force)`; `force = -dU/dr`. No unit conversion done here —
callers pass already-converted params. `r` may be scalar or `np.ndarray`.

- `exp(param_list, r)` — `U = A·exp(-alpha·r)`; params `[A, alpha, …]`.
- `srd(param_list, r)` — Short-Range Damped: `U = P1/(r^|P2| + P3^|P2|)`. Guards
  divide-by-zero in the force.
- `shtr(param_list, r)` — Shifted-Truncated: defined for `r ≤ P3`, shifted so
  `U(P3)=0` and `U'(P3)=0`. Uses `temp_r = min(r, P3)`.
- `powe(param_list, r)` — power law `U = P1·r^P2`; params `[P1, P2, …]`.
- `quarbond(param_list, r)` — quartic bond
  `U = (P2/2)(r-P1)² + (P3/3)(r-P1)³ + (P4/4)(r-P1)⁴`; params `[P1,P2,P3,P4]`.

## Empty-structure builders
- `gen_empty_bonded()` — returns the canonical empty bonded dict for one molecule:
  `ATO{All,Virtual}`, `BON{HAR,QUA}`, `ANG{HAR,QUA}`, `BD3{QBB,MUB}`,
  `DIH{HAR,NCO,COS}`, `CDI{CNCO,CCOS}`, `EXC[]`. (Used by parsing **and** by
  `populate_bonded.build_new_molecule_bonded`.)
- `gen_empty_nonbonded()` — returns `{InteractionType: []}` for all nonbonded types.
- `_params_dict(dic, *args)` — ensures `dic[args]` exists (empty list); used to
  key bonded terms by their parameter tuple. Mutates in place.

## `.off` section parsing
- `_find_off_keywords(off_file_str)` — splits the whole file into the 5 named
  sections by locating their header strings. Returns `{section: substring}`.
- `_recognize_keywords(section)` — regex-finds all `[KEYWORD]` tokens and returns
  `[[keyword, start, end], …]` (matches 3- and 4-letter keywords like `CDIH`).
- `_filter_interactions(interactions)` — splits the keyword list into
  `(bonded_list, nonbonded_list)` by predefined keyword sets.
- `_find_end_bonded(final_bonded, ff_input)` — returns an `['END', loc, loc]`
  marker at the next `[…]` after the last bonded term (delimits the bonded block).
- `_find_molnames(bonded, ff_input)` — extracts molecule names from each `MOL`
  keyword. Returns `list[str]`.
- `_split_into_molecules(bonded)` — groups the flat bonded keyword list into one
  sublist per molecule (split on `MOL`, `END` markers inserted).
- `_gather_fitted_bonded(ff_input)` — splits the intra-potential section into one
  string per fitted bonded interaction block.

## Bonded interaction parsing
- `_parse_bonded(unsorted_bonded, molecule, ff_input)` — orchestrates parsing of a
  single molecule; starts from `gen_empty_bonded()` and fills each section via
  `_parse_bonded_section`. Returns the populated bonded dict.
- `_parse_bonded_section(uns_bonded, section, beginning, ending, ff_input)` —
  the heavy lifter; per `section` keyword parses lines + pulls fitted params from
  `uns_bonded` (advancing `total_bonded_added`). Handles:
  - `ATO` → `{All, Virtual}` (virtual atoms have `*`).
  - `BON`/`ANG` → `{HAR, QUA}` keyed by param tuple → list of atom-index lists.
  - `BD3` → `{QBB, MUB}` (marked *needs testing*).
  - `DIH` → `{HAR, NCO, COS}`.
  - `CDIH` → `{CNCO, CCOS}` (marked *needs testing*).
  - `EXC` → list of int lists.

## Nonbonded parsing + filtering
- `_clean_inter_potential(inter_potential)` — parses the inter-potential section
  into `[(At1,At2)_sorted, 'InteractionType', P1, P2, …]` with floats; strips colons.
- `_filter_bonded(bonded)` — removes empty interaction terms from a molecule's
  bonded dict; passes `EXC` through; drops `Virtual` if empty. Deep-copies.
- `_remove_empty_and_cou_interactions_nonbonded(nonbonded)` — drops empty entries
  and all `COU` interactions; returns a new dict.
- `_remove_netf_torq_atname(all_atoms)` — list of atom **names** excluding NETF/TORQ.
- `_remove_netf_torq_atnum(all_atoms)` — list of atom **numbers** (int) excluding NETF/TORQ.
