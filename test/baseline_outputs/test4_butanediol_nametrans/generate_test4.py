#!/usr/bin/env python3
"""
Test 4: Butanediol - Name Translation Test
Tests name_translation dictionary to map atom names from .off file to topology file
"""
import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))

import afmtogmx as afm

# Change to test directory
test_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(test_dir)

# Define name translation dictionary
# Maps .off file atom names -> topology file atom names
name_translation = {
    'C1': 'CA',
    'C2': 'CB',
    'O1': 'OH',
    'H1': 'HO',
    'H2': 'HA',
    'H3': 'HB'
}

# Read .off file
off = afm.ReadOFF(off_loc="../../sample_off_files/butanediol_intra.off")

# Generate nonbonded tabulated potentials
nonbonded_tabpot = off.gen_nonbonded_tabpot()
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot,
                            prefix='table',
                            name_translation=name_translation)

# Generate bonded tabulated potentials
bonded_tabpot = off.gen_bonded_tabpot()
off.write_bonded_tabpot(bonded_tabpot=bonded_tabpot, prefix='table')

# Generate topology files with name translation
off.gen_nonbonded_topology(template_file='template.top',
                           write_to='temp_nonbonded.top',
                           name_translation=name_translation)
off.gen_bonded_topology(template_file='temp_nonbonded.top',
                        write_to='topol.top',
                        name_translation=name_translation,
                        bonded_tabpot=bonded_tabpot)

print("Test 4 baseline generation complete!")
print(f"Name translation: {name_translation}")
