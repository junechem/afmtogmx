#!/usr/bin/env python3
"""
Test 3: Water - Manual Charge Assignment Test
Tests that manually assigned charges appear correctly in the topology file
Focuses on H20QM molecule type only
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
off = afm.ReadOFF(off_loc="../../sample_off_files/water_intra.off")

# Manually assign charges to the H20QM molecule
# TIP3P-like charges: O=-0.82, H=+0.41, E=0.0 (neutral)
if 'H20QM' in off.charges:
    off.charges['H20QM']['OQM'] = -0.82
    off.charges['H20QM']['HQM'] = 0.41
    off.charges['H20QM']['EQM'] = 0.0

# Generate nonbonded tabulated potentials (only H20QM molecule)
nonbonded_tabpot = off.gen_nonbonded_tabpot(incl_mol=['H20QM'])
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, prefix='table')

# Generate bonded tabulated potentials (only H20QM molecule)
bonded_tabpot = off.gen_bonded_tabpot(incl_mol=['H20QM'])
off.write_bonded_tabpot(bonded_tabpot=bonded_tabpot, prefix='table')

# Generate topology files
off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top', incl_mol=['H20QM'])
off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top',
                        incl_mol=['H20QM'], bonded_tabpot=bonded_tabpot)

print("Test 3 baseline generation complete!")
print(f"Assigned charges: {off.charges.get('H20QM', {})}")
