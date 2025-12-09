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
    -   [ ] Document `exp` function.
    -   [ ] Document `srd` function.
    -   [ ] Document `shtr` function.
    -   [ ] Document `powe` function.
    -   [ ] Document `quarbond` function.
    -   [ ] Document `gen_empty_bonded` function.
    -   [ ] Document `gen_empty_nonbonded` function.
    -   [ ] Document `_find_off_keywords` function.
    -   [ ] Document `_recognize_keywords` function.
    -   [ ] Document `_filter_interactions` function.
    -   [ ] Document `_find_end_bonded` function.
    -   [ ] Document `_find_molnames` function.
    -   [ ] Document `_split_into_molecules` function.
    -   [ ] Document `_gather_fitted_bonded` function.
    -   [ ] Document `_parse_bonded` function.
    -   [ ] Document `_parse_bonded_section` function.
    -   [ ] Document `_clean_inter_potential` function.
    -   [ ] Document `_filter_bonded` function.
    -   [ ] Document `_remove_empty_and_cou_interactions_nonbonded` function.
    -   [ ] Document `_remove_netf_torq_atname` function.
    -   [ ] Document `_remove_netf_torq_atnum` function.

3.  **`src/afmtogmx/core/tabulated_potentials.py` - Potential Calculations**
    -   [ ] Document key public/complex private functions.

4.  **`src/afmtogmx/core/topology.py` - Topology File Generation**
    -   [ ] Document key public/complex private functions.

5.  **`src/afmtogmx/core/residues.py` - Residue Logic**
    -   [ ] Document key functions.

6.  **`src/afmtogmx/core/compare.py` - Comparison Logic**
    -   [ ] Document key functions.

7.  **Final Review**
    -   [ ] Review all new docstrings for consistency and clarity.
