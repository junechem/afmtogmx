#!/usr/bin/env python3
"""
Test 2: Ethane - Bonded Interactions Test
Tests workflow with bonded interactions (bonds, angles, dihedrals)
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
off = afm.ReadOFF(off_loc="../../sample_off_files/ethane_intra.off")

# Generate nonbonded tabulated potentials
nonbonded_tabpot = off.gen_nonbonded_tabpot()
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, prefix='table')

# Generate bonded tabulated potentials
bonded_tabpot = off.gen_bonded_tabpot()
off.write_bonded_tabpot(bonded_tabpot=bonded_tabpot, prefix='table')

# Generate topology files
off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')
off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')

print("Test 2 baseline generation complete!")
