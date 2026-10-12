# `report.py` — publication-style text export

Renders the parsed force-field data on a `ReadOFF` as a fixed-width plain-text
document, sized for pasting into Word at default margins with Courier New 10pt.
**Units are the `.off` file's native kcal / Ångström — nothing is converted.**

Public entry point is `write_report(off, path, ...)`; `ReadOFF.write_report`
(`gen_md.py`) is a thin delegate that returns `self`. Everything else in the
module is a private renderer.

## Module constants

- `LINE_WIDTH = 95` — banner width and the target maximum line length.
- `ATOM_W = 6` — column width for atom-name cells.
- `BONDED_SCHEMAS` — `{category: {subtype: schema}}` where `category` is
  `BON`/`ANG`/`BD3`/`DIH`/`CDI` and a schema is
  `{'natoms', 'formula', 'cols': [(header, fmt_spec), ...], 'remap'?}`.
  `remap` (present on `DIH/NCO` and `DIH/COS`) reorders the stored tuple before
  formatting: stored `(V_D, m, delta)` → rendered `(delta, V_D, m)`.
  `BONDED_SCHEMAS['CDIH']` is aliased to `BONDED_SCHEMAS['CDI']` because
  `_parse_bonded` stores under `CDIH` while `gen_empty_bonded` preallocates `CDI`.
- `BONDED_CATEGORY_LABEL` — category key → display label ("Bonds", "Angles", …).
- `BONDED_CATEGORY_ORDER = ['BON', 'ANG', 'BD3', 'DIH', 'CDI', 'CDIH']` — render order.
- `NONBONDED_ORDER` — section render order, Coulombics first then Table-6 vdW
  forms: `COU, THC, GLJ, BUC, DBU, STR, EXP, GEX, POW, CSP, GDP, PEX, DPO, SRD, TTP`.
  Prefixes present in the data but absent from this list are appended in
  first-appearance order.

## Public entry point

`write_report(off, path, incl_mol=None, name_translation=None, molname_translation=None, notation='standard')`
- Renders every molecule in `off.bonded` that passes `incl_mol` (falsy `incl_mol`
  ⇒ all), joins the lines, and writes them to `path` (overwritten).
- `name_translation` = `{off_atom_name: display_name}`, applied to every atom cell.
  `molname_translation` = `{off_molecule_name: display_name}`, applied to banners
  and section headers. Unlisted names render unchanged.
- **Validates while rendering.** Unknown bonded subtypes and unresolvable atom
  numbers are collected into an `issues` list; if it is non-empty at the end,
  raises `ValueError` ("write_report refusing to emit: …") and **writes nothing**.
- `notation` must be `'standard'` or `'PN'`, else `ValueError` — checked before
  any work is done.

## The two notations

Selected by `notation`, and applied only to the **nonbonded** sections (bonded
column labels come from `BONDED_SCHEMAS` and never change).

- `'standard'` — publication-style labels from the manual (`A`, `alpha`, `Cn`,
  `R0`, `beta`, `r0`, …), one hand-written layout per 3-char prefix.
- `'PN'` — generic `P1, P2, …` sized to the longest row, no value transformation.

For `POW`, `SRD`, `STR`, `PEX`, `GDP`, `DPO`, `CSP` and `TTP` the `'standard'`
layout applies `abs()` to the stored exponent so the `n` column shows a
magnitude; **the coefficient keeps its fitted sign**, so `U = Cn / r^n` with
positive `n` still recovers the correct attractive/repulsive contribution. The
coefficient column header is the literal string `Cn(kcal/mol)` — the exponent is
not substituted into the label (see commit `90ea8c9`).

## Layout resolution

- `_prefix(itype)` — first 3 characters of a stored interaction-type key, which
  is what names the functional form per the CRYOFF manual (Table 6 footnote,
  p. 23). Trailing characters are the fitter's own identifier, so `EXPINTRA`,
  `EXPINTER` and `EXPW` all render under `EXP`, and `STRC` under `STR`.
  Non-string keys are returned unchanged.
- `_nonbonded_layout(prefix, rows, notation)` → `(formula, [(header, fmt), …], row_transform)`.
  One `if` branch per supported prefix. `row_transform` is the identity for most
  forms; the length-guarded `abs()` transforms described above for the rest.
  `notation='PN'` short-circuits to `_pn_layout`, as does an unrecognized prefix.
- `_pn_layout(prefix, rows)` — generic `P1..Pn` columns sized to the longest row.
- `_fmt_width(spec)` — pulls the width digits out of a format spec (`'>16.3f'` → `16`)
  so headers align with their values.

## Per-molecule rendering

- `_render_molecule(off, molname, name_tr, mol_tr, issues, notation)` — assembles
  one molecule's block: banner → charges → bonded → nonbonded → closing rule.
  Emits `(no bonded interactions)` when nothing rendered.
- `_render_charges(charges_for_mol, name_tr)` — `Type / Partial-charges(e)` table;
  returns `[]` when the molecule has no charge entries.

## Bonded rendering

- `_render_bonded_category(category_dict, category, ato_map, molname, name_tr, issues)`
  — loops the category's subtypes, skipping empty ones. A subtype with no schema
  is **not** rendered; it appends an issue instead, which later aborts the write.
- `_render_bonded_subtype(...)` — header + column line + rows. Atom numbers are
  resolved to names via `ato_map`, then **deduplicated** by
  `(canonical_signature, param_tuple)`, so symmetry-equivalent orderings of the
  same term appear once.
- `_canonical_atom_signature(names)` — the dedupe key. 2 atoms: order-independent.
  3 atoms: central atom fixed, ends sorted. 4 atoms: equal to its reverse.
  Other lengths: the tuple as-is.
- `_atom_name(ato_map, atnum)` — resolves an atom number through `ATO['All']`,
  trying both `int` and `str` keys (the parser path determines which), and
  returns the **Coulomb type** (`entry[1]`). Raises `KeyError(atnum)` if absent.
- `category_of(label)` — inverse of `BONDED_CATEGORY_LABEL`, used only to build
  error messages. Public by name but not part of the intended API.

## Nonbonded rendering

- `_select_nonbonded_pairs(nonbonded, mol_atom_names)` — keeps a pair if **either**
  atom belongs to the molecule, so solute–solvent cross terms survive `incl_mol`
  filtering and appear in the included molecule's section.
- `_render_nonbonded(pairs, name_tr, notation)` — discovers which prefixes are
  present, orders them by `NONBONDED_ORDER` (extras appended in first-appearance
  order), resolves each layout **once**, prints an `Energy expressions:` block of
  the formulas, then renders each section. Returns
  `['(no nonbonded interactions)', '']` when `pairs` is empty.
- `_render_nonbonded_section(prefix, pairs, name_tr, layout)` — one row per stored
  parameter set whose 3-char prefix matches. The formula is omitted here (it is
  already in the `Energy expressions:` block) to keep lines within `LINE_WIDTH`.
  Rows shorter than the column list are padded with blanks rather than raising.
  Emits `(none for selected molecule)` if the section matched nothing.
