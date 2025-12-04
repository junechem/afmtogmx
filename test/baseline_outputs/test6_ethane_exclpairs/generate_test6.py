#!/usr/bin/env python3
"""
Test 6: Ethane - Excluded Pairs Test
Tests excl_pairs parameter to exclude specific atom pair interactions
"""
import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))

import afmtogmx as afm

# Change to test directory
test_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(test_dir)

# Define excluded pairs (exclude C-C and H-H interactions)
excl_pairs = [['C', 'C'], ['H', 'H']]

# Read .off file
off = afm.ReadOFF(off_loc="../../sample_off_files/ethane_intra.off")

# Generate nonbonded tabulated potentials with excluded pairs
nonbonded_tabpot = off.gen_nonbonded_tabpot(excl_pairs=excl_pairs)
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, prefix='table')

# Generate bonded tabulated potentials
bonded_tabpot = off.gen_bonded_tabpot()
off.write_bonded_tabpot(bonded_tabpot=bonded_tabpot, prefix='table')

# Generate topology files with excluded pairs
off.gen_nonbonded_topology(template_file='template.top',
                           write_to='temp_nonbonded.top',
                           excl_pairs=excl_pairs)
off.gen_bonded_topology(template_file='temp_nonbonded.top',
                        write_to='topol.top',
                        bonded_tabpot=bonded_tabpot)

print("Test 6 baseline generation complete!")
print(f"Excluded pairs: {excl_pairs}")
