#!/usr/bin/env python3
"""
Test 8: Ethane - Clean Configuration-Based Workflow
Demonstrates how the configuration system enables cleaner, more maintainable code
"""
import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))

import afmtogmx as afm

# Change to test directory
test_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(test_dir)

print("=" * 70)
print("Test 8: Clean Configuration-Based Workflow")
print("=" * 70)

# Read .off file and configure all defaults in one place
print("\nStep 1: Initialize and configure")
print("-" * 70)
off = afm.ReadOFF(off_loc="../../sample_off_files/ethane_intra.off")

# Set all default parameters once
off.set_config(
    spacing_nonbonded=0.0005,
    length_nonbonded=3.0,
    spacing_bonded=0.0001,
    length_bonded=0.3,
    tabpot_prefix='table',
    tabpot_dir='',
    scale_C6=True
)

print("[OK] ReadOFF initialized and configured")

# Load charges from file using method chaining
print("\nStep 2: Load charges from file")
print("-" * 70)
off.load_charges_from_file('charges.txt')
print("[OK] Charges loaded")
print(f"     C: {off.charges['UNK']['C']}")
print(f"     H: {off.charges['UNK']['H']}")

# Generate nonbonded potentials (uses config values, stores automatically)
print("\nStep 3: Generate nonbonded tabulated potentials")
print("-" * 70)
off.gen_nonbonded_tabpot()
print("[OK] Nonbonded potentials generated and stored")

# Write nonbonded potentials (uses stored potentials)
print("\nStep 4: Write nonbonded tabulated potentials")
print("-" * 70)
off.write_nonbonded_tabpot()
print("[OK] Nonbonded potentials written")

# Generate bonded potentials (uses config values, stores automatically)
print("\nStep 5: Generate bonded tabulated potentials")
print("-" * 70)
off.gen_bonded_tabpot()
print("[OK] Bonded potentials generated and stored")

# Write bonded potentials (uses stored potentials)
print("\nStep 6: Write bonded tabulated potentials")
print("-" * 70)
off.write_bonded_tabpot()
print("[OK] Bonded potentials written")

# Generate topology files
print("\nStep 7: Generate topology files")
print("-" * 70)
off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')
print("[OK] Nonbonded topology generated")

off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')
print("[OK] Bonded topology generated")

print("\n" + "=" * 70)
print("Test 8 Complete!")
print("=" * 70)
print("\nThis test demonstrates the clean configuration-based workflow:")
print("  1. Set all defaults once with set_config()")
print("  2. Generate methods store results automatically")
print("  3. Write methods use stored results - no variables needed!")
print("  4. Use method chaining for cleaner code")
print("\nCompare this to traditional approach (see test2) where you would:")
print("  - Assign results to variables: nonbonded_tabpot = off.gen_nonbonded_tabpot(...)")
print("  - Pass variables back: off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot)")
print("  - Repeat spacing, length, prefix, etc. for every method call")
print("\nGenerated files:")
print("  - topol.top")
print("  - temp_nonbonded.top")
print("  - tabpot/table_*.xvg (bonded and nonbonded)")
print("=" * 70)
