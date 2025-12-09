"""The afmtools_v2 package contains tools related to the production of input files for a gromacs
simulation based on a .off file produced from CRYOFF. The main functions of the program are found
within the package gen_md. We suggest using the package in the following way:

    >>> import afmtogmx as afm
    >>> off = afm.ReadOFF(off_loc = "path/to/intra.off")

From here, bonded and nonbonded parameters are stored as dictionaries at off.bonded, off.nonbonded


CONFIGURATION SYSTEM
====================

The ReadOFF class supports a configuration system that allows you to set default parameters for all workflow methods.
This is especially useful when you want to use the same parameter values across multiple method calls.

To set configuration values:

    >>> off.set_config(spacing_nonbonded=0.001, scale_C6=False, tabpot_prefix='my_table')

To get configuration values:

    >>> all_config = off.get_config()  # Returns full config dictionary
    >>> spacing = off.get_config('spacing_nonbonded')  # Returns specific value

Configuration keys include: incl_mol, excl_interactions, excl_pairs, spacing_nonbonded, length_nonbonded,
scale_C6, sc_sigma, spacing_bonded, length_bonded, tabpot_prefix, tabpot_dir, write_blank, name_translation, special_pairs

When calling workflow methods, explicitly passed parameters override config values, which override built-in defaults.


CHARGE ASSIGNMENT
=================

If charges are needed, they can be manually assigned to the off.charges dictionary. By default, all atom charges are
0.0. The dictionary format is: {"MolName": {"AtomName": charge_value}}. For example:

    >>> off.charges['H20QM']['OQM'] = -0.82
    >>> off.charges['H20QM']['HQM'] = 0.41

Alternatively, charges can be loaded from a file using load_charges_from_file():

    >>> off.load_charges_from_file('charges.txt')

The file format should be:
    MOLNAME1
    Atom1 Charge1
    Atom2 Charge2
    MOLNAME2
    Atom3 Charge3

Any atoms not specified in the file will remain at their default charge of 0.0.


GENERATING TABULATED POTENTIALS
================================

To generate nonbonded tabulated potentials, the recommended workflow is as follows:

    >>> off.gen_nonbonded_tabpot()  # Stores result in off.nonbonded_tabpot

Several options are available for gen_nonbonded_tabpot(), but by default, tabulated potentials are generated for all
nonbonded interactions (besides coulombic) between atoms. The result is automatically stored in off.nonbonded_tabpot
for convenient access. You can still capture the return value if needed for advanced use cases.
See the helpdocs for more in depth discussion


To generate bonded tabulated potentials, the recommended workflow is as follows:

    >>> off.gen_bonded_tabpot()  # Stores result in off.bonded_tabpot

Several more options are available; see the helpdocs for more in depth discussion


To produce topology files, use functions gen_nonbonded_topology() and gen_bonded_topology(), along with a template.off
file. An example for the template.off file can be seen in the afmtools_v2/examples directory located on your device.
Recommended workflow for generating nonbonded topology files is

    >>> off.gen_nonbonded_topology(template_file = 'template.top', write_to = 'temp_nonbonded.top')


Recommended workflow for generating bonded topology files (assuming nonbonded .top file 'temp_nonbonded.top' was
created) is

    >>> off.gen_bonded_topology(template_file = 'temp_nonbonded.top', write_to='topol.top')

Several more options are available; see the helpdocs for more in depth discussion


Finally, to write out tabulated potentials, use the functions write_nonbonded_tabpot() and write_bonded_tabpot().
These methods automatically use the stored potentials from the gen methods:

    >>> off.write_nonbonded_tabpot()  # Uses off.nonbonded_tabpot
    >>> off.write_bonded_tabpot()     # Uses off.bonded_tabpot

You can still pass custom dictionaries explicitly if needed for advanced use cases.

So, a good final workflow going from intra.off file to functional force field is:

    >>> import afmtogmx as afm
    >>> off = afm.ReadOFF(off_loc = "path/to/intra.off")
    >>> # Manually set charges if needed:
    >>> # off.charges['MolName']['AtomName'] = charge_value
    >>> # Or load from file:
    >>> # off.load_charges_from_file('charges.txt')
    >>> off.gen_nonbonded_tabpot()  # Stores in off.nonbonded_tabpot
    >>> off.gen_bonded_tabpot()     # Stores in off.bonded_tabpot (optional)
    >>> off.write_nonbonded_tabpot()  # Uses off.nonbonded_tabpot
    >>> off.write_bonded_tabpot()     # Uses off.bonded_tabpot
    >>> off.gen_nonbonded_topology(template_file = 'template.top', write_to = 'temp_nonbonded.top')
    >>> off.gen_bonded_topology(template_file = 'temp_nonbonded.top', write_to='topol.top')


Recommended: Using the configuration system for cleaner code:

    >>> import afmtogmx as afm
    >>> off = afm.ReadOFF(off_loc = "path/to/intra.off")
    >>> off.set_config(tabpot_prefix='my_table', tabpot_dir='tables/')
    >>> off.load_charges_from_file('charges.txt')
    >>> off.gen_nonbonded_tabpot()  # Stores automatically, uses config values
    >>> off.gen_bonded_tabpot()     # Stores automatically, uses config values
    >>> off.write_nonbonded_tabpot()  # Uses stored result and config prefix/dir
    >>> off.write_bonded_tabpot()     # Uses stored result and config prefix/dir
    >>> off.gen_nonbonded_topology(template_file = 'template.top', write_to = 'temp_nonbonded.top')
    >>> off.gen_bonded_topology(template_file = 'temp_nonbonded.top', write_to='topol.top')

"""
from afmtogmx.core import topology, tabulated_potentials, gen_md, functions, compare
from afmtogmx.core.gen_md import ReadOFF
