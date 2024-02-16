"""The afmtools_v2 package contains tools related to the production of input files for a gromacs
simulation based on a .off file produced from CRYOFF. The main functions of the program are found
within the package gen_md. We suggest using the package in the following way:

    >>> import afmtogmx as afm
    >>> off = afm.ReadOFF(off_loc = "path/to/intra.off"

From here, bonded and nonbonded parameters are stored as dictionaries at off.bonded, off.nonbonded


If charges are available in the .off file, they can be calculated using

    >>> off.calc_charges(Options)

This populates the dictionary stored at off.charges with atom charges. By default, all atoms charges are 0.0 until
off.calc_charges(Options) is called. The dictionary off.charges may be set by hand, but it must match the format of the
off.charges dictionary
See the helpdocs for what options may be used with off.calc_charges(Options))


To generate nonbonded tabulated potentials, the recommended workflow is as follows:

    >>> nonbonded_tabpot = off.gen_nonbonded_tabpot()

Several options are available for gen_nonbonded_tabpot(), but by default, tabulated potentials are generated for all
nonbonded interactions (besides coulombic) between atoms. A dictionary containing tables is returned to allow
flexibility with scaling columns by certain factors, for example for HFE calculations.
See the helpdocs for more in depth discussion


To generate bonded tabulated potentials, the recommended workflow is as follows:

    >>> bonded_tabpot = off.gen_bonded_tabpot()

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

    >>> off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot)
    >>> off.write_bonded_tabpot(bonded_tabpot=bonded_tabpot)

So, a good final workflow going from intra.off file to functional force field is:

    >>> import afmtogmx as afm
    >>> off = afm.ReadOFF(off_loc = "path/to/intra.off")
    >>> off.calc_charges(Options)  # only use this line if charges present in .off file
    >>> nonbonded_tabpot = off.gen_nonbonded_tabpot()
    >>> bonded_tabpot = off.gen_bonded_tabpot()  # optional depending on if bonded tabulated interactions in .off file
    >>> off.write_nonbonded_tabpot(nonbonded_tabpot=nonbonded_tabpot)
    >>> off.write_bonded_tabpot(bonded_tabpot=bonded_tabpot)
    >>> off.gen_nonbonded_topology(template_file = 'template.top', write_to = 'temp_nonbonded.top')
    >>> off.gen_bonded_topology(template_file = 'temp_nonbonded.top', write_to='topol.top', bonded_tabpot = bonded_tabpot)  # bonded_tabpot is optional depending on if there are bonded tabulated interactions defined in the .off file

"""
from afmtogmx.core import topology, tabulated_potentials, gen_md, functions, compare
from afmtogmx.core.gen_md import ReadOFF
