# `gmx_backend.py` — `GROMACSBackend` (`off.gmx`)

The public GROMACS API. Accessed as `off.gmx` on a `ReadOFF`. Owns GROMACS config
and the generated tabulated potentials; reads shared FF data from
`self._parent.{bonded, nonbonded, charges}`. Delegates real work to
`tabulated_potentials.py` and `topology.py`.

## Construction / state
`__init__(self, parent)` — stores the parent `ReadOFF`, initializes `self.config`
(see keys below) and `self.nonbonded_tabpot = None`, `self.bonded_tabpot = None`.

**Config keys / defaults:** `special_pairs={}`, `incl_mol=[]`,
`excl_interactions=[]`, `excl_pairs=[]`, `spacing_nonbonded=0.0005`,
`length_nonbonded=3`, `scale_C6=True`, `sc_sigma=0.0`, `spacing_bonded=0.0001`,
`length_bonded=0.3`, `tabpot_prefix='table'`, `tabpot_dir=''`, `write_blank=True`,
`name_translation={}`.

**Parameter resolution** (used throughout via a `SimpleNamespace` comprehension):
explicit arg (if not `None`) → `self.config` value → built-in default.

## Config methods
- `set_config(**kwargs)` — updates `self.config` in place; returns `self` (chainable).
- `get_config(key=None)` — returns one value, or a copy of the whole config dict.

## Tabulated potential generation
- `gen_nonbonded_tabpot(special_pairs=None, incl_mol=None, excl_interactions=None,
  excl_pairs=None, spacing_nonbonded=None, length_nonbonded=None, scale_C6=None,
  sc_sigma=None)` — builds nonbonded tables. Pipeline:
  `_gen_included_atoms` → `_filter_nonbonded` → `_gen_nonbond_tabpam`; if
  `sc_sigma != 0`, additionally `_scale_for_FE` using a generated nonbonded string.
  Stores + returns the dict in `self.nonbonded_tabpot`.
- `gen_bonded_tabpot(incl_mol=None, spacing_bonded=None, length_bonded=None)` —
  builds bonded tables per molecule via `gen_bonded_tabpam`, threading `num_tables`
  so indices stay unique. Empty `incl_mol` ⇒ all molecules. Stores + returns
  `self.bonded_tabpot` (`{molname: {param_tuple: [table_no, x, U, F]}}`).

## Topology generation
- `gen_nonbonded_topology(name_translation=None, template_file=None, incl_mol=None,
  excl_interactions=None, excl_pairs=None, scale_C6=None, special_pairs=None,
  write_to=None, sc_sigma=None)` — inserts a generated `[ nonbond_params ]` block
  into `template_file` and writes it. `template_file` is **required** (returns with
  a message if missing/unreadable). Default `write_to` is `nonbond_topol.top` beside
  the template. With `sc_sigma != 0`, rewrites each line's C12 to `C6·sc_sigma⁶`.
- `gen_bonded_topology(name_translation=None, incl_mol=None, template_file="",
  write_to="", bonded_tabpot=None)` — generates per-molecule bonded sections and
  splices them into `template_file`'s `[ moleculetype ]` blocks; writes the result.
  Uses `self.bonded_tabpot` if `bonded_tabpot` not passed (needed for QUA/QBB).
  Default `write_to` is `bonded_topol.top` beside the template. Empty `incl_mol` ⇒ all.

## Writing tables to disk
- `write_nonbonded_tabpot(nonbonded_tabpot=None, tabpot_prefix=None, tabpot_dir=None,
  name_translation=None, write_blank=None)` — writes one `.xvg` per pair (falls back
  to `self.nonbonded_tabpot`; messages if none). Prints the unique-atom list (for
  `energygrps`) and pair list (for `energygrp_table`). Writes the blank table if
  `write_blank`.
- `write_bonded_tabpot(bonded_tabpot=None, tabpot_prefix=None, tabpot_dir=None)` —
  writes one `.xvg` per bonded table (falls back to `self.bonded_tabpot`; messages
  if none).
