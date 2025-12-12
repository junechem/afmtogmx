"""
The afmtogmx package contains tools related to the production of input files for a GROMACS
simulation based on a .off file produced from CRYOFF. The main functions of the program are found
within the core modules. We suggest using the package in the following way:

    >>> import afmtogmx
    >>> off = afmtogmx.ReadOFF(off_loc="path/to/intra.off")

From here, bonded and nonbonded parameters are stored as dictionaries at off.bonded, off.nonbonded.

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
scale_C6, sc_sigma, spacing_bonded, length_bonded, tabpot_prefix, tabpot_dir, write_blank, name_translation, special_pairs.

When calling workflow methods, explicitly passed parameters override config values, which override built-in defaults.

CHARGE ASSIGNMENT
=================

By default, all atom charges are 0.0. If charges are needed, they can be assigned to the `off.charges` dictionary:

    >>> off.charges['H20QM']['OQM'] = -0.82
    >>> off.charges['H20QM']['HQM'] = 0.41

Alternatively, charges can be loaded from a file using `load_charges_from_file()`:

    >>> off.load_charges_from_file('charges.txt')

The file format should be:
    MOLNAME1
    Atom1 Charge1
    Atom2 Charge2

Any atoms not specified in the file will remain at their default charge of 0.0.

EXAMPLE WORKFLOW
================

A typical workflow from an `.off` file to a complete GROMACS force field is as follows:

    >>> import afmtogmx
    >>>
    >>> # 1. Initialize the reader with the .off file
    >>> off = afmtogmx.ReadOFF(off_loc="path/to/intra.off")
    >>>
    >>> # 2. Set any desired global configurations
    >>> off.set_config(tabpot_dir='tables/', name_translation={'C_off': 'C_top'})
    >>>
    >>> # 3. Load charges if necessary
    >>> off.load_charges_from_file('charges.txt')
    >>>
    >>> # 4. Generate all required tabulated potentials (stored internally)
    >>> off.gen_nonbonded_tabpot()
    >>> off.gen_bonded_tabpot()
    >>>
    >>> # 5. Write the tabulated potential .xvg files to the configured directory
    >>> off.write_nonbonded_tabpot()
    >>> off.write_bonded_tabpot()
    >>>
    >>> # 6. Generate the nonbonded topology from a template
    >>> off.gen_nonbonded_topology(template_file='template.top', write_to='temp_nonbonded.top')
    >>>
    >>> # 7. Generate the final, complete topology
    >>> off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')

"""

__version__ = "0.1.0"

from .core.gen_md import ReadOFF

__all__ = [
    "ReadOFF",
]
