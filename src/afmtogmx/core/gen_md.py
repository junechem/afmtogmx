from afmtogmx.core import topology, tabulated_potentials, functions, chargefxns
""" This module contains the main class, ReadOFF, which is used to generate input files for gmx
"""


class ReadOFF:
    """This class reads an .off file
    """

    def __init__(self, off_loc):
        """Initialize the ReadOFF class by reading an intra.off file at location off_loc

        Assuming you create an object using
            off = ReadOFF(off_loc = 'intra.off')

        You can access the completed bonded and nonbonded parameters with the dictionaries stored at
            off.bonded
            off.nonbonded

        The charges (which are all zero by default) can be accessed with
            off.charges

        :param off_loc: string, location of intra.off file to read.
        """
        self.off_loc = off_loc
        self._ff_bonded = {}  # Need to read bonding information at top of .off file to understand how to parse fitted
        # bonded parameters
        self.bonded = {}  # Dictionary to contain all fitted bonded information for each molecule
        self.nonbonded = {}
        self.sections = ""  # dictionary with 5 sections: ff_input, intra_potential, inter_potential,
        # molecular_definition, and table_potential

        self._gen_sections_dict()  # Calls funtion to generate sections dict
        self._gen_bonded()  # Creates self.bonded dictionary with all sections populated with fited parameters
        self.charges = chargefxns._gen_empty_charge_dict(self.bonded)  # Charges dictionary with each atom set to 0.0
        # charge. Using calc_charges will update self.charges with the proper charges, or self.charges can be manually
        # set
        self._gen_nonbonded()  # Creates self.nonbonded dictionary with all sections populated with fitted parameters
        self.residues = {k : None for k, v in self.bonded.items()}

    def _gen_sections_dict(self):
        """Loads an off file into memory, breaks into sections, and stores it as the variable self.sections"""
        try:
            off = open(self.off_loc, 'r').read()
        except FileNotFoundError as e:
            print(f"{e}")
            print("Problem in _gen_sections_dict")
            exit(1)

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

    def calc_charges(self, known_atom=None, known_atom_charge=None, normalization="M-POPULOUS", known_charge_sign=None,
                     tolerance=1E-5):
        """Populates self.charges with nonzero charges derived from the .off file. Will always round charges to 5 digits

        :param known_atom: atom which other charges should be derived from
        :param known_atom_charge: if charge of known_atom is known (eg BLYPSP-4F proton charge), calculate charges based on known_atom and known_atom_charge
        :param normalization: method of normalizition/neutralizing molecules. Currently only support M-POPULOUS: see below
        :param known_charge_sign: Either '+' or '-'. Sign of known_atom charge; to be used in conjuction with known_atom iff known_atom_charge is not used.
        :param tolerance: If excess total charge > tolerance, perform normalization

        Methods of Normalization:
        M-POPULOUS:     For each molname, the number of unique atoms per molecule are counted, and the most populous
        atom (with n atoms) with a nonzero charge has the excess-total-charge/n subtracted from the charge, resulting in
        functionally zero charge.
        """
        charges = chargefxns._gen_empty_charge_dict(self.bonded)

        if known_atom and known_atom_charge:
            charges = chargefxns._gen_charges_from_known(charges, self.nonbonded, known_atom, known_atom_charge)
        elif known_atom and known_charge_sign:
            known_atom_charge = chargefxns._get_known_atom_charge(known_atom, self.nonbonded, sign=known_charge_sign)
            charges = chargefxns._gen_charges_from_known(charges, self.nonbonded, known_atom, known_atom_charge)

        self.charges = chargefxns._normalize_charges(normalization=normalization.upper(), charges=charges,
                                                    bonded=self.bonded, tolerance=tolerance)

    def gen_nonbonded_tabpot(self, special_pairs={}, incl_mol=[], excl_interactions=[], spacing=0.0005, length=3,
                             scale_C6=True):
        """Return dictionary holding nonbonded tabulated potentials for all pairs. By default, it is assumed that the
        only attractive interactions are 'POW_6', 'DPO_6', 'SRD_6', 'PEX_6', that is, POW interactions which are raised
        to the 6th power, DPO interactions which are raised to the 6th power, and so on.

        gen_nonbonded_tabpot can be used with several options including special_pairs, incl_mol, excl_interactions, and
        scale_C6.

            * special_pairs is a non-required argument; if used, it should be a dictionary holding information following
             the format: {pair : list}. The list should contain interactions which will be put in columns 4 and 5 in the
             tabulated potential for that pair. All other interactions found in the .off file for that pair will be
             placed in the 6th and 7th (repulsive) columns. For example, if it is desired that a POW_6 and POW_8
             interaction be placed in the 4 and 5 columns of a tabulated potential (as opposed to the default behavior
             of putting POW_8 in the 6 and 7 columns) for pair ('At1', 'At2'), the user can call:

            gen_nonbonded_tabpot(special_pairs = {('At1', 'At2') : ['POW_6', 'POW_8']}, scale_C6 = False)

            Note 'scale_C6 = False'; when using custom attractive interactions, scaling the C6 coefficient to enable
            the use of dispersion corrections is not supported.

            * incl_mol is a non-required argument; if used, it should be a list of molnames for which tabulated
            potentials should be generated for. For example, if an .off file contains 6 molnames, MOL1, MOL2,...,MOL6,
            and tabulated potentials with only MOL1 and MOL2 are desired, this can be achieved with:

            gen_nonbonded_tabpot(incl_mol = ['MOL1', 'MOL2'])


            * excl_interactions is a non-required argument; if used, it should be a list of interactions to not write
             tabulated potentials for. For example, if there is a solvent buckingham interaction BUC (as is the case for
             WAIL, rWAIL, BLYPSP-4F) that one does not want to write tabulated potentials for, this interaction can be
             excluded with:

            gen_nonbonded_tabpot(excl_interactions = ['BUC']

            Note that the interaction must be written exactly as it is written in the .off file, i.e., if the buckingham
             water interaction is specified as BUCWATER, then to exclude this, we should specify

             excl_interactions = ['BUCWATER']


            *  scale_C6 is a required boolean argument; by default it is set to True. This scales columns 4 and 5
             for each pair that has an attractive potential by C6; for all default attractive interactions, this is
             parameter 2 as listed in the CRYOFF manual. scale_C6 is disabled for custom attractive interactions.

        :param special_pairs: Non-Required. dictionary following format {pair : list} where 'list' is a list of
        interactions that should be considered attractive and placed in columns 4 and 5 for that pair
        :param incl_mol: Non-Required. list containing molnames which tabulated potentials should be generated for.
        Default behavior is to write tabulated potentials for all molnames
        :param excl_interactions:  Non-required. list containing interactions EXACTLY as they appear in the .off file
        which tabulated parameters should not be written for.
        :param spacing:  Default: 0.0005 nm. Float (in nm) for spacing between x-values in the generated tables.
        :param length:   Default: 3 nm. Float (in nm) for the total length of the generated tables
        :param scale_C6: Default: True. Boolean which controls whether columns 4 and 5 of generated tables are scaled
         to allow for dispersion corrections in gromacs simulations.
        :return: dictionary containing {pair : [x_values, COU_pot, COU_force, ATT_pot, ATT_force, REP_pot, REP_force]}
        for each pair; the elements in the list for each pair are numpy arrays.
        """
        total_incl_atoms = tabulated_potentials._gen_included_atoms(incl_mol, self.bonded)  # list containing atoms which are in
        # incl_mol
        filtered_nonbonded = tabulated_potentials._filter_nonbonded(nonbonded=self.nonbonded, excluded_int=excl_interactions,
                                                                    incl_atoms=total_incl_atoms)

        return tabulated_potentials._gen_nonbond_tabpam(nonbonded=filtered_nonbonded, spec_pairs=special_pairs, spacing=spacing,
                                                        length=length, scale_C6=scale_C6)

    def gen_bonded_tabpot(self, incl_mol=[], spacing=0.0001, length=0.3):
        """Generates bonded tabulated potentials for each molecule containing bond type 'QUA' in the .off file.

        Returned Dictionary Format: {'Molname' : {(parameters) : [table_number, x_values, U(x) values, F(x) values], ...}

        :param incl_mol: list, Non-required; list containing molnames to exclusively generate bonded tabualted parameters for
        :param spacing: float, Default = 0.0001 nm; Spacing between x-values in final table
        :param length: float/int, Default = 0.3 nm; Total length of tables to generate
        :return: dictionary
        """
        bonded_tabpams = dict()  # Holds final tabpam dictionary
        num_tables = 0  # Keep track of number of tables that will be written
        if not incl_mol:  # If incl_mol not specified, try to generate tables for all molecules
            incl_mol = self.bonded.keys()

        for molname in self.bonded:
            if molname in incl_mol:
                # generate tabulated parameters
                bonded_tabpams[molname], num_tables = tabulated_potentials.gen_bonded_tabpam(bond_dict=self.bonded[molname],
                                                                                             spacing=spacing, length=length,
                                                                                             num_tables=num_tables)

        return bonded_tabpams  # return completed tabpam dictionary

    def gen_nonbonded_topology(self, name_translation={}, template_file="", incl_mol=[], excl_interactions=[],
                               scale_C6=True, special_pairs={}, write_to=""):
        """Generates [ nonbond_params ] section of topology file and writes topology file (default name nonbond_topol.top).

        :param name_translation: Dictionary, optional; format {'AtIn.off' : 'AtIn.top',...}; If atom not in name_translation, atom name from .off will be used
        :param template_file: Path/To/Template.top, required; should contain at least 1 section [ nonbond_params ] with a blank line following it.
        :param incl_mol: List of strings, optional; molnames to write nonbond_params for.
        :param excl_interactions: List of strings, optional; interactions as written in the .off file which should not be included in the nonbond_params section. Use case: Specific BLYPSP-4F interactions, interactions that should be written to [ pairs ], etc.
        :param scale_C6: Boolean, optional, default = True; adjust C6 parameters to allow for correct dispersion corrections in the gromacs simulation
        :param special_pairs: Dictionary, optional; see discussion in gen_nonbonded_tabpot
        :param write_to: Path/to/output.top, optional; Default behavior is to write to nonbond_topol.top in current directory
        """
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

        incl_atoms = tabulated_potentials._gen_included_atoms(incl_mol, self.bonded)  # generate atoms to include
        filtered_nonbonded = tabulated_potentials._filter_nonbonded(nonbonded=self.nonbonded, excluded_int=excl_interactions,
                                                                    incl_atoms=incl_atoms)  # filter based on excl, incl
        template_nonbonded_location = topology._find_keyword_location(template_file=f, keyword="nonbond_params",
                                                                      begin=0)  # find location of [nonbonded_params ]
        #  in template, taking acount of comments, etc.
        nonbonded_string = topology._gen_nonbonded_string(scale_C6=scale_C6, special_pairs=special_pairs,
                                                          name_translation=name_translation,
                                                          nonbonded=filtered_nonbonded)  # generate string containing
        #  all nonbonded pairs with proper C6, C12
        with open(write_to, 'w') as w:  # Write to file location
            w.write(f[:template_nonbonded_location[1]])
            w.write(nonbonded_string)
            w.write(f[template_nonbonded_location[1]:])

    def gen_bonded_topology(self, name_translation={}, incl_mol=[], template_file="", write_to="", bonded_tabpot={}):
        """Writes bonded portion of topology file.

        * name_translation dictionary allows atom names to be converted from the .off name to the desired name in a .top file. Format is {'.offAtom1' : '.topAtom1', ...}

        * incl_mol list will write bonded section only for certain molnames included in the list. Must match molname in .off and .top file

        * template_file string is a path to the template file to be used. File must include [ moleculetype ] with molecule names matching those in the .off file, and each  [ keyword ] section that should be populated with information from the .off file

        * write_to string is a path to the file that should be written with completed bonded sections. Default is current directory, file named 'bonded_topol.top'

        * bonded_tabpot dictionary is required when: 1) Quartic bonds exist in the .off file 2) QBB exists in the .off file. This dictionary is output by gen_bonded_tabpot()

        :param name_translation: dictionary, non-required
        :param incl_mol: list, non-required
        :param template_file: string, required
        :param write_to: string, non-required
        :param bonded_tabpot: dictionary, non-required (Except in cases where quartic bonds or QBB interactions exist)
        :return: Writes topology file with completed bonded sections
        """


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
        if not incl_mol:  # if no incl_mol list provided, use all molnames
            incl_mol = self.bonded.keys()
        filtered_bonded = dict()
        for molname in incl_mol:  # remove empty functions from self.bonded
            filtered_bonded[molname] = functions._filter_bonded(self.bonded[molname])
        output_topology_strings = dict()  # prepare dictionary to hold strings for each section of topology file
        for molname in filtered_bonded:  # generate completed section strings for each [ keyword ] for each [moleculetype] in the .off file
            output_topology_strings[molname] = topology._gen_bonded_section_strings(name_translation,
                                                                                    filtered_bonded[molname],
                                                                                    self.charges, molname,
                                                                                    bonded_tabpot)
        moleculetype_locations = topology._find_moleculetypes(f)
        new_topology = ""
        new_topology += topology._gen_bonded_string(topology=f, locations=moleculetype_locations,
                                                    topology_strings_dict=output_topology_strings)

        topology._write_topology(final_topology=new_topology, write_to=write_to)

    @staticmethod
    def write_nonbonded_tabpot(nonbonded_tabpot = {}, prefix = "MOL", to_dir = "", name_translation = {}, write_blank = True):
        """Write nonbonded tabulated potentials based on dictionary nonbonded_tabpot.

        * nonbonded_tabpot dictionary is generated by gen_nonbonded_tabpot(). Self-made dictionaries also work, but they need to match the format of the dictionary generated by gen_nonbonded_tabpot()

        * prefix string is the prefix to the At1_At2.xvg that you would like your files to be named. Default is MOL_At1_At2.xvg

        * to_dir string is the directory which you would like to write tabulated potentials to. Default behavior is to create a directory "tabpot" and populate that with the tabulated potentials.

        * name_translation dictionary is a dictionary with the format {'.offAtom1': '.topAtom1', ...} to translate atom names

        * write_blank boolean is True by default. This writes a blank nonbonded tabulated potential file with the name prefix.xvg, as is required for gromacs simulations to work

        :param nonbonded_tabpot: dictionary, required
        :param prefix: string, optional
        :param to_dir: string, optional
        :param name_translation: dictionary, optional
        :return: None
        """
        if not nonbonded_tabpot:
            print("Must include nonbonded_tabpot dictionary, as generated by gen_nonbonded_tabpot(). "
                  "Self-Made dictionary may be used, but it must match the format of the dictionary produced by "
                  "gen_nonbonded_tabpot")
            return

        write_to = tabulated_potentials._to_dir(write_to_dir=to_dir)

        pair_interactions = []
        unique_atoms_notranslation = []

        for atom_pair, tabpot in nonbonded_tabpot.items():  # Write nonbonded tabulated potentials
            tabulated_potentials._write_nonbonded_pair_tabpot(atom_pair=atom_pair, tabpot = tabpot, name_translation = name_translation, write_to = write_to, prefix = prefix)
            pair_interactions.append(tabulated_potentials._translate_pairs(atom_pair=atom_pair, name_translation=name_translation))
            if atom_pair[0] not in unique_atoms_notranslation:
                unique_atoms_notranslation.append(atom_pair[0])
            if atom_pair[1] not in unique_atoms_notranslation:
                unique_atoms_notranslation.append(atom_pair[1])

        unique_atoms_translated = []
        for atom in unique_atoms_notranslation:  # translate atoms for writing energygrps
            if atom in name_translation:
                unique_atoms_translated.append(name_translation[atom])
            else:
                unique_atoms_translated.append(atom)

        print("\nUnique atoms for energygrps: \n")
        print(' '.join(unique_atoms_translated))

        print("\nTotal pair interactions for energygrp_table: \n")
        individual_pairs = [' '.join(i) for i in pair_interactions]
        print('  '.join(individual_pairs))

        if write_blank:  # write blank file if requested
            tabulated_potentials._write_blank_nonbonded(prefix=prefix, write_to = write_to)

    @staticmethod
    def write_bonded_tabpot(bonded_tabpot = {}, prefix ="MOL", to_dir =""):
        """Generate bonded tabulated potentials for either QUA bonds or QBB terms in the .off file.

        :param bonded_tabpot: dictionary, output by gen_nonbonded_tabpot()
        :param prefix: string, Default "MOL", prefix to add to name of .xvg files
        :param to_dir: string, Default "tabpot" subdirectory of current directory
        :return: None
        """
        if not bonded_tabpot:
            print("Must include nonbonded_tabpot dictionary, as generated by gen_nonbonded_tabpot(). "
                  "Self-Made dictionary may be used, but it must match the format of the dictionary produced by "
                  "gen_nonbonded_tabpot")
            return

        write_to = tabulated_potentials._to_dir(write_to_dir=to_dir)

        for molname, molname_dict in bonded_tabpot.items():
            for params, table_list in molname_dict.items():
                tabulated_potentials._write_bonded_tabpot(tabpot=table_list, prefix=prefix, write_to=write_to)




