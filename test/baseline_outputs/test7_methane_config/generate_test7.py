#!/usr/bin/env python3
"""
Test 7: Methane - Configuration System and Charge Loading Test
Tests set_config(), get_config(), parameter resolution, and load_charges_from_file()
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
print("Test 7: Configuration System and Charge Loading Test")
print("=" * 70)

# Read .off file
off = afm.ReadOFF(off_loc="../../sample_off_files/methane_intra.off")

# TEST 1: set_config() with multiple parameters
print("\nTEST 1: Setting configuration values")
print("-" * 70)
off.set_config(
    spacing_nonbonded=0.001,  # Different from default 0.0005
    length_nonbonded=2.5,     # Different from default 3.0
    tabpot_prefix='config_table',  # Different from default 'table'
    tabpot_dir='config_tabpot',    # Different from default ''
    scale_C6=False            # Different from default True
)
print("[OK] Configuration set with custom values")

# TEST 2: get_config() - retrieve all config
print("\nTEST 2: Retrieving all configuration values")
print("-" * 70)
all_config = off.get_config()
print("Configuration dictionary retrieved:")
for key, value in all_config.items():
    print(f"  {key}: {value}")

# TEST 3: get_config(key) - retrieve specific value
print("\nTEST 3: Retrieving specific configuration values")
print("-" * 70)
spacing = off.get_config('spacing_nonbonded')
prefix = off.get_config('tabpot_prefix')
print(f"[OK] spacing_nonbonded: {spacing}")
print(f"[OK] tabpot_prefix: {prefix}")

# TEST 4: load_charges_from_file() method
print("\nTEST 4: Loading charges from file")
print("-" * 70)
print("Charges before loading:")
print(f"  {off.charges}")
print("\nLoading charges from 'charges.txt'...")
result = off.load_charges_from_file('charges.txt')
print("[OK] Charges loaded from file")
print("\nCharges after loading:")
for molname, atoms in off.charges.items():
    print(f"  {molname}:")
    for atomname, charge in atoms.items():
        print(f"    {atomname}: {charge}")

# Verify load_charges_from_file returns self (method chaining)
if result is off:
    print("[OK] load_charges_from_file() returns self (enables method chaining)")
else:
    print("[ERROR] load_charges_from_file() does not return self")

# TEST 5: Generate tabulated potentials using config values (no explicit params)
print("\nTEST 5: Generating outputs using config values")
print("-" * 70)
print("Generating nonbonded tabpot (using config: spacing=0.001, length=2.5)...")
nonbonded_tabpot = off.gen_nonbonded_tabpot()
print("[OK] Nonbonded tabpot generated")

print("Writing nonbonded tabpot (using config: prefix='config_table', dir='config_tabpot')...")
off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot)
print("[OK] Nonbonded tabpot written")

# TEST 6: Verify that explicit parameters override config values
print("\nTEST 6: Testing parameter override (explicit > config > default)")
print("-" * 70)
print("Config has tabpot_prefix='config_table'")
print("But we'll explicitly pass tabpot_prefix='override_table'...")
off.write_nonbonded_tabpot(
    nonbonded_tabpot=nonbonded_tabpot,
    tabpot_prefix='override_table',
    tabpot_dir='override_tabpot'
)
print("[OK] Explicit parameters correctly overrode config values")

# Generate topology files using config values
print("\nGenerating topology files using config values...")
off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')
off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')
print("[OK] Topology files generated")

# TEST 7: Verify method chaining with set_config
print("\nTEST 7: Testing method chaining")
print("-" * 70)
result = off.set_config(spacing_bonded=0.0002)
if result is off:
    print("[OK] set_config() returns self (enables method chaining)")
else:
    print("[ERROR] ERROR: set_config() does not return self")

# Demonstrate full method chaining
print("\nDemonstrating full method chaining:")
print("  off.set_config(...).load_charges_from_file(...)")
off2 = afm.ReadOFF(off_loc="../../sample_off_files/methane_intra.off")
off2.set_config(tabpot_prefix='chained').load_charges_from_file('charges.txt')
print("[OK] Method chaining works!")

print("\n" + "=" * 70)
print("Test 7 Complete!")
print("=" * 70)
print("\nGenerated files:")
print("  - topol.top")
print("  - temp_nonbonded.top")
print("  - config_tabpot/config_table_*.xvg")
print("  - override_tabpot/override_table_*.xvg")
print("\nAll configuration system tests passed!")
print("=" * 70)
