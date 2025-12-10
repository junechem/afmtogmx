from afmtogmx.core import topology, tabulated_potentials, functions, residues
from types import SimpleNamespace
""" This module contains the main class, ReadOFF, which is used to generate input files for gmx
"""


class ReadOFF:
    """Orchestrates the conversion of a CRYOFF .off file to GROMACS inputs.

    This class serves as the main entry point for the `afmtogmx` library. It
    reads and parses a `.off` file, storing the molecular topology, force
    field parameters, and interaction data. It then provides methods to
    generate GROMACS topology (`.top`) and tabulated potential (`.xvg`) files.

    Parameters
    ----------
    off_loc : str
        The file path to the `.off` file to be processed.

    Attributes
    ----------
    off_loc : str
        Stores the provided location of the `.off` file.
    bonded : dict
        A nested dictionary containing all parsed bonded interactions,
        organized by molecule name, interaction type (e.g., 'BON', 'ANG'),
        and specific parameters.
    nonbonded : dict
        A nested dictionary containing all parsed nonbonded interactions,
        organized by atom type pairs (e.g., ('C', 'H')) and interaction
        type (e.g., 'POW_6').
    charges : dict
        A nested dictionary holding the atomic charges for each atom in each
        molecule. Initialized to zero by default and can be populated from
        an external file.
    config : dict
        A dictionary holding configuration options for file generation, such
        as table spacing, potential scaling, and naming conventions.
    residues : dict
        A dictionary to store residue definitions and assignments.
    nonbonded_tabpot : dict or None
        Stores the generated nonbonded tabulated potential data after a call
        to `gen_nonbonded_tabpot`.
    bonded_tabpot : dict or None
        Stores the generated bonded tabulated potential data after a call
        to `gen_bonded_tabpot`.

    Examples
    --------
    >>> from afmtogmx.core.gen_md import ReadOFF
    >>> off = ReadOFF('path/to/your/forcefield.off')
    >>> print(off.bonded.keys())
    dict_keys(['MOL1'])

    """

    def __init__(self, off_loc):
        self.off_loc = off_loc
        self._ff_bonded = {}  # Need to read bonding information at top of .off file to understand how to parse fitted
        # bonded parameters
        self.bonded = {}  # Dictionary to contain all fitted bonded information for each molecule
        self.nonbonded = {}
        self.sections = ""  # dictionary with 5 sections: ff_input, intra_potential, inter_potential,
        # molecular_definition, and table_potential

        # Initialize configuration with defaults
        self.config = {
            'special_pairs': {},
            'incl_mol': [],
            'excl_interactions': [],
            'excl_pairs': [],
            'spacing_nonbonded': 0.0005,
            'length_nonbonded': 3,
            'scale_C6': True,
            'sc_sigma': 0.0,
            'spacing_bonded': 0.0001,
            'length_bonded': 0.3,
            'tabpot_prefix': 'table',
            'tabpot_dir': '',
            'write_blank': True,
            'name_translation': {}
        }

        # Initialize tabulated potential storage
        self.nonbonded_tabpot = None
        self.bonded_tabpot = None

        self._gen_sections_dict()  # Calls funtion to generate sections dict
        self._gen_bonded()  # Creates self.bonded dictionary with all sections populated with fited parameters
        # Initialize charges dictionary with all charges set to 0.0
        # Format: {"Mol1" : {'At1' : 0.0, 'At2' : 0.0...}, 'Mol2':...}
        # User can manually set charges via self.charges dictionary
        self.charges = {mol: {pair[1]: 0.0 for key, pair in self.bonded[mol]['ATO']['All'].items()
                              if pair[1] != 'NETF' and pair[1] != 'TORQ'}
                        for mol in self.bonded}
        self._gen_nonbonded()  # Creates self.nonbonded dictionary with all sections populated with fitted parameters
        self.residues = {"Definitions" : {k : {'All' : functions._remove_netf_torq_atname(v['ATO']['All'])} for k, v in self.bonded.items()}, "Residues" : {k : {'All' : [functions._remove_netf_torq_atnum(v['ATO']['All'])]} for k, v in self.bonded.items()}}

    def _gen_sections_dict(self):
        """Loads an off file into memory, breaks into sections, and stores it as the variable self.sections"""
        try:
            off = open(self.off_loc, 'r').read()
        except FileNotFoundError as e:
            print(f"{e}")
            print("Problem in _gen_sections_dict")
            raise

        self.sections = functions._find_off_keywords(off_file_str=off)

    def _gen_bonded(self):
        keywords_and_locations = functions._recognize_keywords(section=self.sections['ff_input'])  # find all keywords
        # and locations in the self.sections['ff_input'] part of off file
        bonded, nonbonded = functions._filter_interactions(keywords_and_locations)  # splitting interactions into
        # bonded, nonbonded
        bonded.append(functions._find_end_bonded(bonded[-1], self.sections['ff_input']))  # adding 'END' to bonded
        # section to help with parsing
        molnames = functions._find_molnames(bonded, self.sections['ff_input'])  # list of molnames
        molecules = functions._split_into_molecules(bonded)  # list of list of lists; each element in the list
        # corresponds to a different molecule
        unsorted_fitted_bonded = functions._gather_fitted_bonded(self.sections['intra_potential'])
        for molname, molecule in zip(molnames, molecules):
            self.bonded[f'{molname}'] = functions._parse_bonded(unsorted_fitted_bonded, molecule,
                                                                ff_input=self.sections['ff_input'])
        functions.total_bonded_added = 0

    def _gen_nonbonded(self):  # populates self.nonbonded dictionary with correct pairs and parameters
        cleaned_inter_potential = functions._clean_inter_potential(self.sections['inter_potential'])

        for interaction in cleaned_inter_potential:  # populate actual self.nonbonded dictionary
            atom_pair = interaction[0]
            inter_term = interaction[1]
            if 'COU' in inter_term:
                inter_term = 'COU'
            params = interaction[2:]

            if atom_pair not in self.nonbonded:  # if there is not already an interaction for the atom pair
                self.nonbonded[atom_pair] = dict()  # create empty dict for atom pair
                #                self.nonbonded[atom_pair] = functions.gen_empty_nonbonded()  # create empty dict for atom pair
                self.nonbonded[atom_pair][f'{inter_term}'] = []
                self.nonbonded[atom_pair][f'{inter_term}'].append(params)  # populate with parameters
            else:  # if atom pair is already in self.nonbonded, just add parameters
                if inter_term not in self.nonbonded[atom_pair]:
                    self.nonbonded[atom_pair][f'{inter_term}'] = []
                    self.nonbonded[atom_pair][f'{inter_term}'].append(params)  # populate with parameters
                else:
                    self.nonbonded[atom_pair][f'{inter_term}'].append(params)  # populate with parameters


    def gen_nonbonded_tabpot(self, special_pairs=None, incl_mol=None, excl_interactions=None, excl_pairs=None,
                             spacing_nonbonded=None, length_nonbonded=None, scale_C6=None, sc_sigma=None):
        """Generate nonbonded tabulated potentials for all relevant atom pairs.

        This method calculates and stores tabulated potential energy and force
        curves for nonbonded interactions. The results are stored in the
        `self.nonbonded_tabpot` attribute. By default, interactions assumed to
        be attractive are 'POW_6', 'DPO_6', 'SRD_6', and 'PEX_6'. Parameters
        not explicitly provided to this method will be resolved from the
        `self.config` dictionary.

        Parameters
        ----------
        special_pairs : dict, optional
            A dictionary defining custom attractive interactions. Format:
            `{('At1', 'At2') : ['INTERACTION_TYPE1', 'INTERACTION_TYPE2']}`.
            The specified interaction types will be placed in columns 4 and 5
            of the tabulated potential for that pair. All other interactions
            for that pair will be placed in the 6th and 7th (repulsive) columns.
            Defaults to `{}` (from `self.config`).
        incl_mol : list of str, optional
            A list of molecule names (as found in the .off file) for which
            tabulated potentials should be generated. If `None`, potentials
            are generated for all molecules. Defaults to `[]` (from `self.config`).
        excl_interactions : list of str, optional
            A list of interaction names (exactly as written in the .off file)
            to be excluded from the tabulated potentials. Useful for excluding
            specific intra-molecular or solvent interactions.
            Defaults to `[]` (from `self.config`).
        excl_pairs : list of list of str, optional
            A list of atom pairs (e.g., `[['C', 'H']]`) for which nonbonded
            tabulated potentials should *not* be produced.
            Defaults to `[]` (from `self.config`).
        spacing_nonbonded : float, optional
            The spacing (in nm) between x-values in the generated tables.
            Defaults to `0.0005` nm (from `self.config`).
        length_nonbonded : float, optional
            The total length (in nm) of the generated tables.
            Defaults to `3.0` nm (from `self.config`).
        scale_C6 : bool, optional
            If `True`, columns 4 and 5 for each pair that has an attractive
            potential are scaled by C6. This is essential for enabling
            dispersion corrections in GROMACS.
            Defaults to `True` (from `self.config`).
        sc_sigma : float, optional
            If non-zero, the tabulated potentials will be scaled by the
            appropriate amount to enable free energy calculations using the
            sc-sigma value in the `grompp.mdp` file. Defaults to `0.0`
            (from `self.config`).

        Returns
        -------
        dict
            A dictionary containing the generated nonbonded tabulated
            potentials. The keys are atom pairs (e.g., `('C', 'H')`), and
            the values are lists of NumPy arrays:
            `[x_values, COU_pot, COU_force, ATT_pot, ATT_force, REP_pot, REP_force]`.
            This dictionary is also stored in `self.nonbonded_tabpot`.

        Notes
        -----
        - **Parameter Resolution**: Parameters are resolved in the following order
          of precedence: explicit argument > `self.config` setting > hardcoded default.
        - **`special_pairs`**: When using `special_pairs`, `scale_C6` should generally be set to `False`
          unless you explicitly manage the scaling within your custom attractive interactions.
        - **`excl_interactions`**: Ensure interaction names match exactly as they appear
          in the `.off` file (e.g., 'BUCWATER' instead of 'BUC' if that is how it's named).

        Examples
        --------
        >>> off = ReadOFF('path/to/forcefield.off')
        >>> # Generate nonbonded potentials using default config values
        >>> off.gen_nonbonded_tabpot()
        >>> # Access the generated tables
        >>> print(off.nonbonded_tabpot.keys())

        >>> # Generate with custom settings for a specific pair and no C6 scaling
        >>> off.set_config(scale_C6=False)
        >>> off.gen_nonbonded_tabpot(
        ...     special_pairs={('At1', 'At2'): ['POW_6', 'POW_8']},
        ...     spacing_nonbonded=0.001
        ... )
        """
        # Resolve parameters: explicit value → config → default
        p = SimpleNamespace(**{
            k: v if v is not None else self.config[k]
            for k, v in locals().items() if k in self.config
        })

        total_incl_atoms = tabulated_potentials._gen_included_atoms(p.incl_mol, self.bonded)  # list containing atoms which are in
        # incl_mol
        filtered_nonbonded = tabulated_potentials._filter_nonbonded(nonbonded=self.nonbonded, excluded_int=p.excl_interactions,
                                                                    incl_atoms=total_incl_atoms, excl_pairs=p.excl_pairs)

        tabpot = tabulated_potentials._gen_nonbond_tabpam(nonbonded=filtered_nonbonded, spec_pairs=p.special_pairs,
                                                        spacing=p.spacing_nonbonded, length=p.length_nonbonded,
                                                        scale_C6=p.scale_C6)
        if p.sc_sigma != 0.0:
            nonbonded_string = topology._gen_nonbonded_string(scale_C6=p.scale_C6, special_pairs=p.special_pairs,
                                                              nonbonded=filtered_nonbonded, name_translation={})
            tabpot = tabulated_potentials._scale_for_FE(sc_sigma=p.sc_sigma, nonbonded_string=nonbonded_string, tabpot=tabpot)
        else:
            pass

        # Store result in attribute for convenient access
        self.nonbonded_tabpot = tabpot
        return tabpot


    def gen_bonded_tabpot(self, incl_mol=None, spacing_bonded=None, length_bonded=None):
        """Generate bonded tabulated potentials for specific molecules.

        This method calculates and stores tabulated potential energy and force
        curves for bonded interactions, specifically for 'QUA' bonds and 'QBB'
        terms found in the .off file. The results are stored in the
        `self.bonded_tabpot` attribute. Parameters not explicitly provided
        to this method will be resolved from the `self.config` dictionary.

        Parameters
        ----------
        incl_mol : list of str, optional
            A list of molecule names (as found in the .off file) for which
            bonded tabulated potentials should be generated. If `None`,
            potentials are generated for all molecules.
            Defaults to all molecules in `self.bonded` (from `self.config`).
        spacing_bonded : float, optional
            The spacing (in nm) between x-values in the generated tables.
            Defaults to `0.0001` nm (from `self.config`).
        length_bonded : float or int, optional
            The total length (in nm) of the generated tables.
            Defaults to `0.3` nm (from `self.config`).

        Returns
        -------
        dict
            A dictionary containing the generated bonded tabulated potentials.
            The format is `{'Molname' : {(parameters) : [table_number, x_values, U(x) values, F(x) values], ...}}`.
            This dictionary is also stored in `self.bonded_tabpot`.

        Notes
        -----
        - **Parameter Resolution**: Parameters are resolved in the following order
          of precedence: explicit argument > `self.config` setting > hardcoded default.
        - **Supported Interactions**: Currently supports 'QUA' bonds and 'QBB' terms.

        Examples
        --------
        >>> off = ReadOFF('path/to/forcefield.off')
        >>> # Generate bonded potentials using default config values
        >>> off.gen_bonded_tabpot()
        >>> # Access the generated tables
        >>> print(off.bonded_tabpot.keys())

        >>> # Generate for a specific molecule with custom spacing
        >>> off.gen_bonded_tabpot(incl_mol=['ETHANE'], spacing_bonded=0.0005)
        """
        # Resolve parameters: explicit value → config → default
        p = SimpleNamespace(**{
            k: v if v is not None else self.config[k]
            for k, v in locals().items() if k in self.config
        })

        bonded_tabpams = dict()  # Holds final tabpam dictionary
        num_tables = 0  # Keep track of number of tables that will be written
        if not p.incl_mol:  # If incl_mol not specified, try to generate tables for all molecules
            p.incl_mol = self.bonded.keys()

        for molname in self.bonded:
            if molname in p.incl_mol:
                # generate tabulated parameters
                bonded_tabpams[molname], num_tables = tabulated_potentials.gen_bonded_tabpam(bond_dict=self.bonded[molname],
                                                                                             spacing=p.spacing_bonded, length=p.length_bonded,
                                                                                             num_tables=num_tables)

        # Store result in attribute for convenient access
        self.bonded_tabpot = bonded_tabpams
        return bonded_tabpams  # return completed tabpam dictionary

    def gen_nonbonded_topology(self, name_translation=None, template_file=None, incl_mol=None, excl_interactions=None, excl_pairs=None,
                               scale_C6=None, special_pairs=None, write_to=None, sc_sigma=None):
        """Generates the `[ nonbond_params ]` section and writes the nonbonded topology file.

        This method constructs the `[ nonbond_params ]` section of a GROMACS
        topology file based on the parsed `.off` file data. It then writes
        this section into a provided template file and saves the result to
        a new topology file. Parameters not explicitly provided will be
        resolved from the `self.config` dictionary.

        Parameters
        ----------
        name_translation : dict, optional
            A dictionary for translating atom names from the `.off` file to
            the desired names in the `.top` file. Format:
            `{'.offAtom1' : '.topAtom1', ...}`. If an atom is not in
            `name_translation`, its name from the `.off` file will be used.
            Defaults to `{}` (from `self.config`).
        template_file : str, required
            Path to a GROMACS topology template file. This file must contain
            at least one `[ nonbond_params ]` section header with a blank line
            following it. The generated nonbond parameters will be inserted
            at this location.
        incl_mol : list of str, optional
            A list of molecule names (as found in the .off file) for which
            nonbonded parameters should be included.
            Defaults to `[]` (from `self.config`).
        excl_interactions : list of str, optional
            A list of interaction names (exactly as written in the .off file)
            to be excluded from the `[ nonbond_params ]` section. Useful for
            specific interactions that might be handled in other sections (e.g.,
            `[ pairs ]`). Defaults to `[]` (from `self.config`).
        excl_pairs : list of list of str, optional
            A list of atom pairs (e.g., `[['C', 'H']]`) for which nonbonded
            parameters should *not* be written in the topology file.
            Defaults to `[]` (from `self.config`).
        scale_C6 : bool, optional
            If `True`, C6 parameters will be adjusted to allow for correct
            dispersion corrections in GROMACS simulations.
            Defaults to `True` (from `self.config`).
        special_pairs : dict, optional
            See discussion in `gen_nonbonded_tabpot`. Defaults to `{}`
            (from `self.config`).
        write_to : str, optional
            The path to the output topology file. If `None`, the default is
            `nonbond_topol.top` in the same directory as the `template_file`.
        sc_sigma : float, optional
            If non-zero, the `C12` parameters in the `[ nonbond_params ]`
            section will be scaled by the appropriate amount to enable free
            energy calculations using the sc-sigma value in the `grompp.mdp`
            file. Defaults to `0.0` (from `self.config`).

        Returns
        -------
        None

        Notes
        -----
        - **Parameter Resolution**: Parameters are resolved in the following order
          of precedence: explicit argument > `self.config` setting > hardcoded default.
        - **`template_file`**: The template file is crucial for defining the
          overall structure of the GROMACS topology. Ensure it has the correct
          `[ nonbond_params ]` header.
        - **Output Location**: If `write_to` is not specified, the output file
          `nonbond_topol.top` will be created in the same directory as `template_file`.

        Examples
        --------
        >>> off = ReadOFF('path/to/forcefield.off')
        >>> # Generate a nonbonded topology file using default settings
        >>> off.gen_nonbonded_topology(template_file='template.top')

        >>> # Generate with custom name translation and a specific output file
        >>> name_trans = {'C_off': 'C_top', 'H_off': 'H_top'}
        >>> off.set_config(name_translation=name_trans, scale_C6=False)
        >>> off.gen_nonbonded_topology(
        ...     template_file='template.top',
        ...     write_to='my_nonbonded.top'
        ... )
        """
        # Resolve parameters: explicit value → config → default
        p = SimpleNamespace(**{
            k: v if v is not None else self.config[k]
            for k, v in locals().items() if k in self.config
        })

        if not template_file:  # if no template file provided, exit out of function
            print("Template required for generating nonbonded section of topology file. Please supply a template file")
            return
        else:
            try:
                f = open(template_file, 'r').read()
            except FileNotFoundError as e:  # If can't find template file, exit
                print(e)
                return

        if not write_to:  # If no write location given, use directory of template file and name nonbond_topol.top
            import re
            from os.path import dirname, abspath
            CWD = re.sub(r"\\", "/", dirname(abspath(template_file)))
            write_to = f'{CWD}/nonbond_topol.top'

        incl_atoms = tabulated_potentials._gen_included_atoms(p.incl_mol, self.bonded)  # generate atoms to include
        filtered_nonbonded = tabulated_potentials._filter_nonbonded(nonbonded=self.nonbonded, excluded_int=p.excl_interactions,
                                                                    incl_atoms=incl_atoms, excl_pairs=p.excl_pairs)  # filter based on excl, incl
        template_nonbonded_location = topology._find_keyword_location(template_file=f, keyword="nonbond_params",
                                                                      begin=0)  # find location of [nonbonded_params ]
        #  in template, taking acount of comments, etc.
        nonbonded_string = topology._gen_nonbonded_string(scale_C6=p.scale_C6, special_pairs=p.special_pairs,
                                                          name_translation=p.name_translation,
                                                          nonbonded=filtered_nonbonded)  # generate string containing
        #  all nonbonded pairs with proper C6, C12

        if p.sc_sigma != 0.0:
            new_nonbonded_string = ""
            for line in nonbonded_string.split('\n'):
                entries = line.split()
                if len(entries) == 0:
                    continue
                C6, C12 = float(entries[-2]), float(entries[-1])
                if C6 != 0 and C12 != 0:
                    atom_pair = tuple(entries[:2])
                    new_c12 = C6*(p.sc_sigma**6)
                    new_nonbonded_string += topology.single_nonbonded_pair_string(pair = tuple(atom_pair), name_translation={}, C6=C6, C12=new_c12)
                else:
                    new_nonbonded_string += line + "\n"
            nonbonded_string = new_nonbonded_string

        with open(write_to, 'w') as w:  # Write to file location
            w.write(f[:template_nonbonded_location[1]])
            w.write(nonbonded_string)
            w.write(f[template_nonbonded_location[1]:])

    def gen_bonded_topology(self, name_translation=None, incl_mol=None, template_file="", write_to="", bonded_tabpot=None):
        """Generates bonded sections and writes the bonded topology file.

        This method constructs the bonded sections of a GROMACS topology file
        (e.g., `[ atoms ]`, `[ bonds ]`, `[ angles ]`, `[ dihedrals ]`)
        based on the parsed `.off` file data. It then inserts these sections
        into a provided template file and saves the result to a new topology file.
        Parameters not explicitly provided will be resolved from the
        `self.config` dictionary.

        Parameters
        ----------
        name_translation : dict, optional
            A dictionary for translating atom names from the `.off` file to
            the desired names in the `.top` file. Format:
            `{'.offAtom1' : '.topAtom1', ...}`.
            Defaults to `{}` (from `self.config`).
        incl_mol : list of str, optional
            A list of molecule names (as found in the .off file) for which
            bonded sections should be included. If `None`, sections are
            generated for all molecules.
            Defaults to all molecules in `self.bonded` (from `self.config`).
        template_file : str, required
            Path to a GROMACS topology template file. This file must include
            `[ moleculetype ]` sections with molecule names matching those in
            the `.off` file, and each `[ keyword ]` section (e.g., `[ atoms ]`)
            that should be populated with information from the `.off` file.
        write_to : str, optional
            The path to the output topology file. If `None`, the default is
            `bonded_topol.top` in the same directory as the `template_file`.
        bonded_tabpot : dict, optional
            A dictionary of bonded tabulated potentials, typically generated
            by `gen_bonded_tabpot()`. If `None`, the method attempts to use
            the potentials stored in `self.bonded_tabpot`. This is **required**
            if the `.off` file contains `QUA` bonds or `QBB` interactions, as
            these require references to external tabulated potential files.
            Defaults to `self.bonded_tabpot` if available, otherwise `{}`.

        Returns
        -------
        None

        Notes
        -----
        - **Parameter Resolution**: Parameters are resolved in the following order
          of precedence: explicit argument > `self.config` setting > hardcoded default.
        - **`template_file`**: The template file is crucial for defining the
          overall structure of the GROMACS topology. Ensure it has the correct
          `[ moleculetype ]` and `[ keyword ]` headers.
        - **`bonded_tabpot` Requirement**: If `QUA` bonds or `QBB` interactions are
          present, `bonded_tabpot` must be provided, usually as the output
          of `gen_bonded_tabpot()`.

        Examples
        --------
        >>> off = ReadOFF('path/to/forcefield.off')
        >>> # First generate necessary bonded tabulated potentials
        >>> off.gen_bonded_tabpot()  # Stores in self.bonded_tabpot
        >>> # Then generate the bonded topology
        >>> off.gen_bonded_topology(template_file='temp_nonbonded.top', write_to='topol.top')
        # Uses stored self.bonded_tabpot automatically

        >>> # With name translation and specific molecules
        >>> off.set_config(name_translation={'C_off': 'C_top'})
        >>> off.gen_bonded_topology(
        ...     incl_mol=['ETHANE'],
        ...     template_file='template.top',
        ...     write_to='ethane_bonded.top'
        ... )
        # Still uses self.bonded_tabpot
        """
        # Use stored attribute if no parameter provided
        if bonded_tabpot is None:
            bonded_tabpot = self.bonded_tabpot if self.bonded_tabpot is not None else {}

        # Resolve parameters: explicit value → config → default
        p = SimpleNamespace(**{
            k: v if v is not None else self.config[k]
            for k, v in locals().items() if k in self.config
        })

        if not template_file:  # if no template file provided, exit out of function
            print("Template required for generating nonbonded section of topology file. Please supply a template file")
            return
        else:
            try:
                f = open(template_file, 'r').read()
            except FileNotFoundError as e:  # If can't find template file, exit
                print(e)
                return
        if not write_to:  # If no write location given, use directory of template file and name nonbond_topol.top
            import re
            from os.path import dirname, abspath
            CWD = re.sub(r"\\", "/", dirname(abspath(template_file)))
            write_to = f'{CWD}/bonded_topol.top'
        if not p.incl_mol:  # if no incl_mol list provided, use all molnames
            p.incl_mol = self.bonded.keys()
        filtered_bonded = dict()
        for molname in p.incl_mol:  # remove empty functions from self.bonded
            filtered_bonded[molname] = functions._filter_bonded(self.bonded[molname])
        output_topology_strings = dict()  # prepare dictionary to hold strings for each section of topology file
        for molname in filtered_bonded:  # generate completed section strings for each [ keyword ] for each [moleculetype] in the .off file
            output_topology_strings[molname] = topology._gen_bonded_section_strings(p.name_translation,
                                                                                    filtered_bonded[molname],
                                                                                    self.charges, molname,
                                                                                    bonded_tabpot)
        moleculetype_locations = topology._find_moleculetypes(f)
        new_topology = ""
        new_topology += topology._gen_bonded_string(topology=f, locations=moleculetype_locations,
                                                    topology_strings_dict=output_topology_strings)

        topology._write_topology(final_topology=new_topology, write_to=write_to)

    def write_nonbonded_tabpot(self, nonbonded_tabpot=None, tabpot_prefix=None, tabpot_dir=None, name_translation=None, write_blank=None):
        """Write nonbonded tabulated potentials to `.xvg` files.

        This method takes the generated nonbonded tabulated potentials (either
        from `self.nonbonded_tabpot` or an explicitly provided dictionary) and
        writes them to `.xvg` files in a specified directory. Each atom pair
        will have its own `.xvg` file.

        Parameters
        ----------
        nonbonded_tabpot : dict, optional
            A dictionary of nonbonded tabulated potentials, typically generated
            by `gen_nonbonded_tabpot()`. If `None`, the method attempts to use
            the potentials stored in `self.nonbonded_tabpot`.
        tabpot_prefix : str, optional
            A prefix for the generated `.xvg` filenames (e.g., `'table'` will
            result in files like `table_At1_At2.xvg`). Defaults to `'table'`
            (from `self.config`).
        tabpot_dir : str, optional
            The path to the directory where the `.xvg` files should be written.
            If the directory does not exist, it will be created. If `None` or
            empty, defaults to a `tabpot` subdirectory in the current working
            directory (from `self.config`).
        name_translation : dict, optional
            A dictionary for translating atom names from the `.off` file to
            the desired names in the `.xvg` filenames. Format:
            `{'.offAtom1': '.topAtom1', ...}`. Defaults to `{}` (from `self.config`).
        write_blank : bool, optional
            If `True`, a blank nonbonded tabulated potential file named
            `<tabpot_prefix>.xvg` will also be written. This is often required
            for GROMACS simulations. Defaults to `True` (from `self.config`).

        Returns
        -------
        None

        Notes
        -----
        - **`nonbonded_tabpot` format**: If an external `nonbonded_tabpot` dictionary
          is provided, it must match the format generated by `gen_nonbonded_tabpot()`:
          `{pair : [x_values, COU_pot, COU_force, ATT_pot, ATT_force, REP_pot, REP_force]}`.
        - **Atom names for `energygrps`**: The method prints a list of unique
          translated atom names and pair interactions to the console, which can
          be useful for defining `energygrps` or `energygrp_table` in GROMACS
          `mdp` files.

        Examples
        --------
        >>> off = ReadOFF('path/to/forcefield.off')
        >>> off.gen_nonbonded_tabpot() # Generate potentials, store in self.nonbonded_tabpot
        >>> off.write_nonbonded_tabpot(tabpot_dir='./my_tabpots')
        # Writes files like './my_tabpots/table_C_H.xvg'

        >>> # Using an external dictionary and a custom prefix
        >>> custom_potentials = off.gen_nonbonded_tabpot(scale_C6=False)
        >>> off.write_nonbonded_tabpot(nonbonded_tabpot=custom_potentials, tabpot_prefix='custom_table')
        # Writes files like './tabpot/custom_table_C_H.xvg' (if tabpot_dir is default)
        """
        # Use stored attribute if no parameter provided
        if nonbonded_tabpot is None:
            nonbonded_tabpot = self.nonbonded_tabpot

        # Check if we have potentials to write
        if nonbonded_tabpot is None:
            print("No nonbonded potentials to write. Call gen_nonbonded_tabpot() first or pass nonbonded_tabpot parameter.")
            return

        # Resolve parameters: explicit value → config → default
        p = SimpleNamespace(**{
            k: v if v is not None else self.config[k]
            for k, v in locals().items() if k in self.config
        })

        write_to = tabulated_potentials._to_dir(write_to_dir=p.tabpot_dir)

        pair_interactions = []
        unique_atoms_notranslation = []

        for atom_pair, tabpot in nonbonded_tabpot.items():  # Write nonbonded tabulated potentials
            tabulated_potentials._write_nonbonded_pair_tabpot(atom_pair=atom_pair, tabpot=tabpot, name_translation=p.name_translation, write_to=write_to, prefix=p.tabpot_prefix)
            pair_interactions.append(tabulated_potentials._translate_pairs(atom_pair=atom_pair, name_translation=p.name_translation))
            if atom_pair[0] not in unique_atoms_notranslation:
                unique_atoms_notranslation.append(atom_pair[0])
            if atom_pair[1] not in unique_atoms_notranslation:
                unique_atoms_notranslation.append(atom_pair[1])

        unique_atoms_translated = []
        for atom in unique_atoms_notranslation:  # translate atoms for writing energygrps
            if atom in p.name_translation:
                unique_atoms_translated.append(p.name_translation[atom])
            else:
                unique_atoms_translated.append(atom)

        print("\nUnique atoms for energygrps: \n")
        print(' '.join(unique_atoms_translated))

        print("\nTotal pair interactions for energygrp_table: \n")
        individual_pairs = [' '.join(i) for i in pair_interactions]
        print('  '.join(individual_pairs))

        if p.write_blank:  # write blank file if requested
            tabulated_potentials._write_blank_nonbonded(prefix=p.tabpot_prefix, write_to=write_to)

    def write_bonded_tabpot(self, bonded_tabpot=None, tabpot_prefix=None, tabpot_dir=None):
        """Write bonded tabulated potentials to `.xvg` files.

        This method takes the generated bonded tabulated potentials (either
        from `self.bonded_tabpot` or an explicitly provided dictionary) and
        writes them to `.xvg` files in a specified directory. Each unique
        bonded potential type will have its own `.xvg` file.

        Parameters
        ----------
        bonded_tabpot : dict, optional
            A dictionary of bonded tabulated potentials, typically generated
            by `gen_bonded_tabpot()`. If `None`, the method attempts to use
            the potentials stored in `self.bonded_tabpot`.
        tabpot_prefix : str, optional
            A prefix for the generated `.xvg` filenames (e.g., `'table'` will
            result in files like `table_b0.xvg`). Defaults to `'table'`
            (from `self.config`).
        tabpot_dir : str, optional
            The path to the directory where the `.xvg` files should be written.
            If the directory does not exist, it will be created. If `None` or
            empty, defaults to a `tabpot` subdirectory in the current working
            directory (from `self.config`).

        Returns
        -------
        None

        Notes
        -----
        - **`bonded_tabpot` format**: If an external `bonded_tabpot` dictionary
          is provided, it must match the format generated by `gen_bonded_tabpot()`.

        Examples
        --------
        >>> off = ReadOFF('path/to/forcefield.off')
        >>> # Generate potentials, store in self.bonded_tabpot
        >>> off.gen_bonded_tabpot()
        >>> off.write_bonded_tabpot(tabpot_dir='./my_bonded_tabpots')
        # Writes files like './my_bonded_tabpots/table_b0.xvg'

        >>> # Using an external dictionary and a custom prefix
        >>> custom_bonded_potentials = off.gen_bonded_tabpot(incl_mol=['MOL2'])
        >>> off.write_bonded_tabpot(bonded_tabpot=custom_bonded_potentials, tabpot_prefix='custom_bonded')
        """
        # Use stored attribute if no parameter provided
        if bonded_tabpot is None:
            bonded_tabpot = self.bonded_tabpot

        # Check if we have potentials to write
        if bonded_tabpot is None:
            print("No bonded potentials to write. Call gen_bonded_tabpot() first or pass bonded_tabpot parameter.")
            return

        # Resolve parameters: explicit value → config → default
        p = SimpleNamespace(**{
            k: v if v is not None else self.config[k]
            for k, v in locals().items() if k in self.config
        })

        write_to = tabulated_potentials._to_dir(write_to_dir=p.tabpot_dir)

        for molname, molname_dict in bonded_tabpot.items():
            for params, table_list in molname_dict.items():
                tabulated_potentials._write_bonded_tabpot(tabpot=table_list, prefix=p.tabpot_prefix, write_to=write_to)


    def gen_residues(self, residue_definition={}, residue_atnums={}):
        """Populate `self.residues` with custom residue information.

        This method allows users to define custom residue groupings based on
        atom types and atom numbers within each molecule. It performs checks
        to ensure that the provided residue definitions and atom numbers
        correspond to atoms present in the `.off` file.

        Parameters
        ----------
        residue_definition : dict, optional
            A dictionary defining residues based on atom types. Default is `{}`.
        residue_atnums : dict, optional
            A dictionary defining residues based on atom numbers. Default is `{}`.

        Returns
        -------
        None

        Notes
        -----
        - **`residue_definition` format**:
          `{molname : {Residue1Name : [AtType1, AtType2, ...], Residue2Name : ...}}`
          `molname` is a molecule name from the .off file. `ResidueNName` is the
          desired residue name. The list `[AtType1, AtType2...]` contains strings
          matching atom types found in the residue. You must specify multiple
          copies of each atom type if they are found in your residue.

          Example for Ethane where 'CH3' is a residue:
          `{'Ethane' : {'CH3' : ['C', 'H', 'H', 'H']}}`

        - **`residue_atnums` format**:
          `{molname : {Residue1Name : [[Atnum1, Atnum2, ...], [AtnumA, AtnumB, ...]]}}`
          `molname` is a molecule name from the .off file. `ResidueNName` is the
          desired residue name. The inner list contains information regarding
          atom numbers that are part of the desired residue.

          Example for Ethane 'CH3' residue, given specific atom numbers:
          If the .off file defines atoms as:
          1   C, 2   H, 3   H, 4   H, 5   C, 6   H, 7   H, 8   H, 9   NETF, 10  TORQ
          Then `residue_atnums` could be:
          `{'Ethane' : {'CH3' : [[1, 2, 3, 4], [5, 6, 7, 8]]}}`

        - The method prints messages indicating the start and completion of
          residue generation. Warnings or errors will be printed if molecule
          names, atom types, or atom numbers are not found in the `.off` file.

        Examples
        --------
        >>> off = ReadOFF('path/to/forcefield.off')
        >>> # Define a simple residue for a molecule named 'METHANE'
        >>> res_def = {'METHANE': {'C1': ['C', 'H', 'H', 'H', 'H']}}
        >>> res_atn = {'METHANE': {'C1': [[1, 2, 3, 4, 5]]}}
        >>> off.gen_residues(residue_definition=res_def, residue_atnums=res_atn)
        Generating Residues
        Done generating residues
        """

        print("Generating Residues")

        # Check if molecule names, atom types, and atom numbers provided are actually in the .off file
        residues._check_residue_definitions(bonded = self.bonded, residue_definition = residue_definition)
        residues._check_residue_atnums(bonded = self.bonded, residue_atnums= residue_atnums)

        self.residues = residues._set_residue_definitions(self.residues, residue_definition)
        self.residues = residues._set_residue_atnums(self.residues, residue_atnums)

        print("Done generating residues")

    def set_config(self, **kwargs):
        """Update the internal configuration with new default values.

        This method provides a convenient way to set or override default
        parameters that will be used in subsequent method calls, such as
        `gen_nonbonded_tabpot`. It modifies the `self.config` dictionary
        in-place.

        The method returns the instance of the class (`self`) to allow for
        method chaining.

        Parameters
        ----------
        **kwargs
            Arbitrary keyword arguments corresponding to the configuration
            keys to be updated.

        Returns
        -------
        ReadOFF
            The instance of the class, allowing for method chaining.

        Notes
        -----
        Available configuration keys and their default values:
        - `special_pairs` (dict): `{}`
        - `incl_mol` (list): `[]`
        - `excl_interactions` (list): `[]`
        - `excl_pairs` (list): `[]`
        - `spacing_nonbonded` (float): `0.0005`
        - `length_nonbonded` (int): `3`
        - `scale_C6` (bool): `True`
        - `sc_sigma` (float): `0.0`
        - `spacing_bonded` (float): `0.0001`
        - `length_bonded` (float): `0.3`
        - `tabpot_prefix` (str): `'table'`
        - `tabpot_dir` (str): `''`
        - `write_blank` (bool): `True`
        - `name_translation` (dict): `{}`

        Examples
        --------
        >>> off = ReadOFF('path/to/forcefield.off')
        >>> off.set_config(spacing_nonbonded=0.001, scale_C6=False)
        >>> # The configuration changes are applied to subsequent method calls.
        >>> off.gen_nonbonded_tabpot()
        >>> # To access the generated tables: off.nonbonded_tabpot

        """
        self.config.update(kwargs)
        return self

    def get_config(self, key=None):
        """Retrieve configuration value(s).

        This method allows access to the current configuration settings.
        If a specific key is provided, it returns the value associated with
        that key. If no key is provided, a copy of the entire configuration
        dictionary is returned.

        Parameters
        ----------
        key : str, optional
            The name of the configuration setting to retrieve.
            If `None`, the entire configuration dictionary is returned.

        Returns
        -------
        value or dict
            The value of the specified configuration setting, or a copy of
            the full configuration dictionary if no key is provided.

        Examples
        --------
        >>> off = ReadOFF('path/to/forcefield.off')
        >>> off.get_config('spacing_nonbonded')
        0.0005
        >>> all_config = off.get_config()
        >>> 'scale_C6' in all_config
        True
        """
        if key is None:
            return self.config.copy()
        return self.config.get(key)

    def load_charges_from_file(self, file_path):
        """Load atomic charges from a file into `self.charges`.

        This method reads charges from a specified file and populates the
        `self.charges` dictionary. Charges are assigned based on molecule name
        and atom name. Any atoms not listed in the file will retain their
        default charge of 0.0.

        Parameters
        ----------
        file_path : str
            The path to the file containing the atomic charges.

        Returns
        -------
        ReadOFF
            The instance of the class, allowing for method chaining.

        Warnings
        --------
        This method will overwrite any previously set charges for atoms
        specified in the input file.

        Notes
        -----
        The charge file should follow this format:

        MOLNAME1
        Atom1 Charge1
        Atom2 Charge2
        ...
        MOLNAME2
        Atom3 Charge3
        ...

        Lines starting with '#' or empty lines are ignored.
        If a molecule name from the file is not found in the force field,
        it will be skipped with a warning.
        If an atom name within a molecule is not found, it will be skipped
        with a warning.

        Examples
        --------
        Assuming 'charges.txt' contains:
        UNK
        C 0.1
        H 0.05

        >>> off = ReadOFF('path/to/forcefield.off')
        >>> off.load_charges_from_file('charges.txt')
        >>> print(off.charges['UNK']['C'])
        0.1
        """
        try:
            with open(file_path, 'r') as f:
                current_mol = None
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):  # Skip empty lines and comments
                        continue

                    # Check if this is a molname header (single word) or atom-charge pair (two words)
                    parts = line.split()
                    if len(parts) == 1:
                        # This is a molname header
                        current_mol = parts[0]
                        if current_mol not in self.charges:
                            print(f"Warning: Molecule '{current_mol}' from charge file not found in force field. Skipping.")
                            current_mol = None
                    elif len(parts) == 2:
                        # This is an atom-charge pair
                        if current_mol is None:
                            print(f"Warning: Atom-charge pair '{line}' found before any molecule name. Skipping.")
                            continue

                        atomname, charge = parts[0], float(parts[1])

                        if atomname not in self.charges[current_mol]:
                            print(f"Warning: Atom '{atomname}' not found in molecule '{current_mol}'. Skipping.")
                            continue

                        self.charges[current_mol][atomname] = charge
                    else:
                        print(f"Warning: Unrecognized line format: '{line}'. Skipping.")

        except FileNotFoundError:
            print(f"Error: Charge file '{file_path}' not found.")
            raise
        except ValueError as e:
            print(f"Error parsing charge file: {e}")
            raise

        return self




