# `tabulated_potentials.py` — GROMACS tabulated potentials (`.xvg`)

Generates and writes tabulated potential tables for both nonbonded and bonded
interactions. Called by `GROMACSBackend` (`off.gmx`). Units are converted here
(`.off` kcal/Å → GROMACS kJ/nm).

## Atom / interaction filtering
- `_gen_included_atoms(incl_mol, bonded)` — returns the list of atom names to
  include. If `incl_mol` given, only those molecules' atoms; else all
  (NETF/TORQ always excluded).
- `_filter_nonbonded(nonbonded, excluded_int, incl_atoms, excl_pairs)` — removes
  excluded interactions/pairs and atoms not in `incl_atoms`; renames variable-power
  interactions to `TYPE_POWER` (e.g. `POW` → `POW_6`). Returns the cleaned dict
  ready for table generation. (`POW/PEX/DPO/SRD` are the variable-power set.)

## Nonbonded table generation
- `_gen_nonbond_tabpam(nonbonded, spec_pairs, spacing, length, scale_C6, charges={})`
  — core generator. For each pair, sums attractive vs repulsive contributions over
  `x = arange(0, length+spacing, spacing)` (x[0] cloned to x[1] to dodge r=0
  infinities). Returns `{pair: [x, COU_pot, COU_force, ATT_pot, ATT_force, REP_pot, REP_force]}`.
  - Default attractive interactions: `POW_6, DPO_6, SRD_6` (columns 4-5).
  - `BUC` is split: attractive `POW^-6` part + repulsive `EXP` part.
  - `THC` is **not supported** (prints a message).
  - `spec_pairs` overrides which interactions are attractive for a pair.
  - **Exits** if `spec_pairs` used with `scale_C6=True`, or if >1 attractive
    interaction occurs with `scale_C6=True`.
- `gen_nonbond_table(x, interaction, params, ATT, scale_C6)` — dispatcher for a
  single interaction. Applies unit conversion then calls the matching
  `functions.*`: `EXP→exp`, `STR→shtr`, `POW→powe`, `SRD→srd`. When `ATT and
  scale_C6`, the leading coefficient is forced to -1 (C6 scaling handled in `.top`).
  `PEX`/`DPO` are **not implemented** (return `(1, 2)`); unknown types print a warning.
- `_scale_for_FE(sc_sigma, nonbonded_string, tabpot)` — scales the repulsive
  columns (`tabpot[pair][5]`/`[6]`) by `1/(C6·sc_sigma⁶)` for free-energy runs that
  use `sc-sigma` in the mdp. Reads C6/C12 from the `[ nonbond_params ]` string.

## Bonded table generation
- `gen_bonded_tabpam(bond_dict, spacing, length, num_tables)` — generates tables
  for `BON['QUA']` quartic bonds and `BD3['QBB']` terms via `functions.quarbond`.
  Returns `({param_tuple: [table_number, x, U, F]}, updated_num_tables)`.
  Returns `({}, num_tables)` if neither QUA nor QBB present. `num_tables` threads
  across molecules to keep table indices unique.

## Output (file writers)
- `_to_dir(write_to_dir="")` — ensures the output directory exists (creates it, or
  a default `./tabpot/`). Returns the path.
- `_write_nonbonded_pair_tabpot(atom_pair, tabpot, name_translation, write_to, prefix)`
  — writes one pair's 7-column table to `<write_to>/<prefix>_<At1>_<At2>.xvg`
  (atom names translated via `name_translation`).
- `_write_blank_nonbonded(prefix, write_to)` — writes the all-zeros default table
  GROMACS requires, `<prefix>.xvg`, out to 5.0 nm @ 0.0005 nm spacing.
- `_write_bonded_tabpot(tabpot, prefix, write_to)` — writes one bonded table
  `[x, U, F]` to `<write_to>/<prefix>_b<table_number>.xvg`.
- `_translate_pairs(atom_pair, name_translation)` — returns the pair with names
  translated (original kept if not in the map). Used to print energygrp tables.
