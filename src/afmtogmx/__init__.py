"""
Convert CRYOFF `.off` force fields into GROMACS and OpenMM input files.

Parsing a `.off` file gives you a `ReadOFF` object holding the force field data,
plus one backend per simulation engine:

    >>> import afmtogmx
    >>> off = afmtogmx.ReadOFF(off_loc="path/to/intra.off")
    >>>
    >>> off.bonded, off.nonbonded, off.charges, off.residues   # shared parsed data
    >>> off.gmx        # GROMACSBackend — tabulated potentials (.xvg) + topology (.top)
    >>> off.openmm     # OpenMMBackend  — force field XML + processed PDB

Methods that modify the force field itself live on `off` and affect both backends;
engine-specific options live on the backend.

CHARGE ASSIGNMENT
=================

All charges default to 0.0. Assign them directly:

    >>> off.charges['H2OQM']['OQM'] = -0.82
    >>> off.charges['H2OQM']['HQM'] = 0.41

or load them from a file (blank lines and '#' comments ignored):

    >>> off.load_charges_from_file('charges.txt')

    MOLNAME1
    Atom1 Charge1
    Atom2 Charge2

An `atom charge` line appearing before any molecule name is applied to every molecule
containing that atom. Atoms not listed keep their default of 0.0.

CONFIGURATION
=============

Each backend keeps its own config, so defaults are set once rather than repeated on
every call. Resolution order is: explicit argument -> config value -> built-in default.

    >>> off.gmx.set_config(spacing_nonbonded=0.001, scale_C6=False, tabpot_dir='tables/')
    >>> off.gmx.get_config()             # full dict
    >>> off.gmx.get_config('scale_C6')   # one value

`off.gmx` keys: spacing_nonbonded, length_nonbonded, spacing_bonded, length_bonded,
scale_C6, sc_sigma, tabpot_prefix, tabpot_dir, write_blank, incl_mol, excl_interactions,
excl_pairs, name_translation, special_pairs.

`off.openmm` keys: incl_mol, molname_translations.

GROMACS WORKFLOW
================

    >>> import afmtogmx
    >>>
    >>> off = afmtogmx.ReadOFF(off_loc="path/to/intra.off")
    >>> off.load_charges_from_file('charges.txt')
    >>> off.gmx.set_config(tabpot_dir='tables/', name_translation={'C_off': 'C_top'})
    >>>
    >>> # generated tables are stored on the backend, so nothing needs passing around
    >>> off.gmx.gen_nonbonded_tabpot()      # -> off.gmx.nonbonded_tabpot
    >>> off.gmx.gen_bonded_tabpot()         # -> off.gmx.bonded_tabpot
    >>> off.gmx.write_nonbonded_tabpot()
    >>> off.gmx.write_bonded_tabpot()
    >>>
    >>> # topology is built in two passes over a template
    >>> off.gmx.gen_nonbonded_topology(template_file='template.top',
    ...                                write_to='temp_nonbonded.top')
    >>> off.gmx.gen_bonded_topology(template_file='temp_nonbonded.top',
    ...                             write_to='topol.top')

OPENMM WORKFLOW
===============

    >>> off = afmtogmx.ReadOFF(off_loc="path/to/intra.off")
    >>> off.load_charges_from_file('charges.txt')
    >>>
    >>> off.openmm.set_config(molname_translations={'H2OQM': 'SOL'})
    >>> off.openmm.gen_xml(output='forcefield.xml')
    >>>
    >>> # rename PDB atoms to match the XML and emit fresh CONECT records
    >>> off.openmm.preprocess_pdb('conf.pdb', 'forcefield.xml', 'conf_processed.pdb')

The XML is self-contained: no second force field file to merge at runtime. OpenMM itself
is an optional dependency, used only to look up atomic masses — without it the XML is
still written, but every `mass` attribute is 0.0.

MODIFYING THE FORCE FIELD
=========================

`change_molecule` replaces one molecule's parameters with a reference force field stored
in `afmtogmx/reference_ff/`. Cross-term (solute-solvent) pairs keep their fitted values:

    >>> off.change_molecule(mol_name='H2OQM', reference_ff='BLYPSP-4F',
    ...                     ref_mol_name='H2OQM')

`populate_bonded` rebuilds a molecule's bonded section from a topology directory of
atoms and bonds, deriving angles, dihedrals and exclusions from the bond graph. It
returns a *new* ReadOFF; the original is untouched:

    >>> new_off = off.populate_bonded(directory='cb7_topology/', new_mol_name='CB7',
    ...                               remove_molecules=['C7F'])

`gen_residues` splits a molecule into named residues by atom type or atom number:

    >>> off.gen_residues(residue_definition={'PROT': {'ALA': ['CA', 'CB', 'HA']}})

`write_report` writes a fixed-width text summary of every fitted parameter in the
`.off` file's native kcal/Angstrom units:

    >>> off.write_report('forcefield_parameters.txt')

See README.md for full examples and docs/ for reference material.
"""

__version__ = "0.1.0"

from .core.gen_md import ReadOFF

__all__ = [
    "ReadOFF",
]
