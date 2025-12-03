#!/usr/bin/env python3
"""
Test 1: Methane - Basic Test
Tests basic workflow with no special options
"""
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))

import afmtogmx as afm

# Change to test directory
test_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(test_dir)

print("=" * 60)
print("Test 1: Methane - Basic Test")
print("=" * 60)

# Read OFF file
off = afm.ReadOFF(off_loc="../../sample_off_files/methane_intra.off")

# Generate nonbonded tabulated potentials
print("\nGenerating nonbonded tabulated potentials...")
nonbonded_tabpot = off.gen_nonbonded_tabpot()

# Write nonbonded tabulated potentials
print("Writing nonbonded tabulated potentials...")
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot, prefix='table')

# Generate nonbonded topology
print("Generating nonbonded topology...")
off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')

# Generate bonded topology
print("Generating bonded topology...")
off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')

print("\n" + "=" * 60)
print("Test 1 Complete!")
print("Generated files:")
print("  - topol.top")
print("  - temp_nonbonded.top")
print("  - tabpot/table_*.xvg")
print("=" * 60)
