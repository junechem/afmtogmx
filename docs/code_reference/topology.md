# `topology.py` — GROMACS `.top` topology generation

Builds the `[ nonbond_params ]` section and the per-molecule bonded sections
(`[ atoms ]`, `[ bonds ]`, `[ angles ]`, `[ dihedrals ]`, etc.) and splices them
into a user-supplied template `.top` file. Called by `GROMACSBackend` (`off.gmx`).
Units converted here (`.off` kcal/Å → GROMACS kJ/nm).

## Template / keyword location
- `_find_keyword_location(template_file, keyword="", begin=0)` — finds the char
  span of a `[ keyword ]` section up to the next blank line, skipping commented
  (`;`) headers and any match before `begin`. **Exits** if not found.
- `_find_moleculetypes(topology)` — finds all uncommented `[ moleculetype ]`
  sections. Returns `[((begin, end), uncommented_text), …]`.

## Nonbonded section
- `_gen_nonbonded_string(scale_C6, special_pairs, name_translation, nonbonded)` —
  computes C6/C12 per pair and returns the full `[ nonbond_params ]` body string.
  - `BUC` → C6 from `|param·4.184·1e-6|`, C12=1 (has repulsive exp part).
  - Repulsive-only / non-default-attractive → C12=1.
  - Default-attractive (`POW_6/DPO_6/SRD_6`) → C6 from `|param·4.184·1e-6|`.
  - `special_pairs` pairs get C6=C12=1.0.
  - **Exits** if `special_pairs` used with `scale_C6=True`, or >1 attractive with
    `scale_C6=True`. `THC` not supported.
- `single_nonbonded_pair_string(pair, name_translation, C6, C12)` — formats one
  `At1 At2 1 C6 C12` line (names translated).

## Bonded sections (per molecule)
- `_gen_bonded_section_strings(name_translation, filtered_bonded, charges, molname, bonded_tabpot)`
  — top-level: builds a skeleton then fills each GROMACS section string from the
  filtered bonded dict. Dispatches ATO→atoms (+ virtual_sitesn), BON→bonds,
  ANG→angles, BD3→bonds+angles, DIH→dihedrals, CDI→cmap, EXC→exclusions.
- `_gen_bonded_skeleton(filtered_bonded)` — returns `{gromacs_keyword: ""}` for
  exactly the sections present (via a fixed `.off`→`.top` keyword map).
- `_gen_bonded_atoms(name_translation, charges, int_dict, molname)` — `[ atoms ]`
  body: atnum, name, resnum(=1), resname(=molname), type, charge group, charge.
  Skips NETF/TORQ.
- `_gen_bonded_bonds(int_dict, molname, bonded_tabpot)` — `[ bonds ]` body.
  `HAR` → func type 1 with converted params; `QUA` → func type 8 referencing the
  table number from `bonded_tabpot[molname][params][0]`. **Exits** with a clear
  message if the QUA table is missing (mismatch between `gen_bonded_tabpot` and
  `gen_bonded_topology`, usually an `incl_mol` inconsistency).
- `_gen_bonded_angles(int_dict)` — `[ angles ]` body; `HAR` func type 1, `QUA`
  func type 6.
- `_gen_bonded_bd3(int_dict, molname="", bonded_tabpot={}, bonds=False)` — BD3
  (3-center) terms. `QBB` writes two table bonds (func 8) when `bonds=True`, else a
  func-3 angle; `MUB` writes a func-4 angle (params flagged uncertain in a print).
- `_gen_dihedrals(int_dict)` — `[ dihedrals ]` body; `HAR`→func 2, `NCO`/`COS`→func 9.
- `_gen_cmap(int_dict)` — `[ cmap ]` body from CDI data.
- `_gen_exclusions(int_dict)` — `[ exclusions ]` body (skips atom 0).

## Assembly + output
- `_gen_bonded_string(topology, locations, topology_strings_dict)` — rebuilds the
  whole `.top`: template header + each `[ moleculetype ]` with its generated
  sections inserted, preserving non-`.off` content.
- `_gen_molname_bonded(topology, bonded)` — inserts a single molecule's generated
  section strings into its `[ moleculetype ]` block at the right keyword headers,
  handling comments and the `system`/`molecules` tail.
- `_write_topology(final_topology, write_to)` — writes the assembled string to file.
