#!/usr/bin/env python3
"""
Test 5: Methane - Soft-Core Sigma Test
Tests sc_sigma parameter for free energy calculations
"""
import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))

import afmtogmx as afm

# Change to test directory
test_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(test_dir)

# Read .off file
off = afm.ReadOFF(off_loc="../../sample_off_files/methane_intra.off")

# Use sc_sigma for soft-core scaling (typical value for free energy calculations)
sc_sigma_value = 0.3

# Generate nonbonded tabulated potentials with sc_sigma
nonbonded_tabpot = off.gen_nonbonded_tabpot(sc_sigma=sc_sigma_value)
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, prefix='table')

# Generate topology files with sc_sigma
off.gen_nonbonded_topology(template_file='template.top',
                           write_to='temp_nonbonded.top',
                           sc_sigma=sc_sigma_value)
off.gen_bonded_topology(template_file='temp_nonbonded.top',
                        write_to='topol.top')

print("Test 5 baseline generation complete!")
print(f"sc_sigma value: {sc_sigma_value}")
