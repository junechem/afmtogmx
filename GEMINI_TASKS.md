# Gemini Tasks: afmtogmx Documentation Overhaul

This file tracks the plan for improving the documentation of the `afmtogmx` codebase, following the NumPy docstring standard.

## Plan

1.  **`src/afmtogmx/core/gen_md.py` - The Core API**
    -   [x] Rewrite `ReadOFF` class docstring.
    -   [x] Rewrite `ReadOFF.__init__` docstring.    -   [x] Rewrite `ReadOFF.set_config` docstring.
    -   [x] Rewrite `ReadOFF.get_config` docstring.
    -   [x] Rewrite `ReadOFF.load_charges_from_file` docstring.
    -   [x] Rewrite `ReadOFF.gen_nonbonded_tabpot` docstring.
    -   [x] Rewrite `ReadOFF.write_nonbonded_tabpot` docstring.
    -   [x] Rewrite `ReadOFF.gen_bonded_tabpot` docstring.
    -   [x] Rewrite `ReadOFF.write_bonded_tabpot` docstring.
    -   [x] Rewrite `ReadOFF.gen_nonbonded_topology` docstring.
    -   [x] Rewrite `ReadOFF.gen_bonded_topology` docstring.
    -   [x] Rewrite `ReadOFF.gen_residues` docstring.

2.  **`src/afmtogmx/core/functions.py` - Low-Level Parsing**
    -   [x] Document `exp` function.
    -   [x] Document `srd` function.
    -   [x] Document `shtr` function.
    -   [x] Document `powe` function.
    -   [x] Document `quarbond` function.
    -   [x] Document `gen_empty_bonded` function.
    -   [x] Document `gen_empty_nonbonded` function.
    -   [x] Document `_find_off_keywords` function.
    -   [x] Document `_recognize_keywords` function.
    -   [x] Document `_filter_interactions` function.
    -   [x] Document `_find_end_bonded` function.
    -   [x] Document `_find_molnames` function.
    -   [x] Document `_split_into_molecules` function.
    -   [x] Document `_gather_fitted_bonded` function.
    -   [x] Document `_parse_bonded` function.
    -   [x] Document `_parse_bonded_section` function.
    -   [x] Document `_clean_inter_potential` function.
    -   [x] Document `_filter_bonded` function.
    -   [x] Document `_remove_empty_and_cou_interactions_nonbonded` function.
    -   [x] Document `_remove_netf_torq_atname` function.
    -   [x] Document `_remove_netf_torq_atnum` function.

3.  **`src/afmtogmx/core/tabulated_potentials.py` - Potential Calculations**
    -   [ ] Document `_gen_included_atoms` function.
    -   [ ] Document `_filter_nonbonded` function.
    -   [ ] Document `_gen_nonbond_tabpam` function.
    -   [ ] Document `gen_nonbond_table` function.
    -   [ ] Document `_scale_for_FE` function.
    -   [ ] Document `gen_bonded_tabpam` function.
    -   [ ] Document `_to_dir` function.
    -   [ ] Document `_write_nonbonded_pair_tabpot` function.
    -   [ ] Document `_write_blank_nonbonded` function.
    -   [ ] Document `_write_bonded_tabpot` function.
    -   [ ] Document `_translate_pairs` function.

4.  **`src/afmtogmx/core/topology.py` - Topology File Generation**
    -   [ ] Document key public/complex private functions.

5.  **`src/afmtogmx/core/residues.py` - Residue Logic**
    -   [ ] Document key functions.

6.  **`src/afmtogmx/core/compare.py` - Comparison Logic**
    -   [ ] Document key functions.

7.  **Final Review**
    -   [ ] Review all new docstrings for consistency and clarity.
