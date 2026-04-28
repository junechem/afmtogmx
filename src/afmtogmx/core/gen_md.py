from afmtogmx.core import functions, residues
from afmtogmx.core.gmx_backend import GROMACSBackend
import warnings
from collections import defaultdict
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
    residues : dict
        A dictionary to store residue definitions and assignments.
    gmx : GROMACSBackend
        GROMACS output backend. Provides all GROMACS-specific methods such
        as ``gen_nonbonded_tabpot`` and topology generation. Configuration
        and generated tabulated potentials are owned by this object.

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
        self.gmx = GROMACSBackend(self)

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

    # ---------------------------------------------------------------------------
    # Deprecated properties — forward to off.gmx equivalents
    # ---------------------------------------------------------------------------

    @property
    def config(self):
        """Deprecated. Use ``off.gmx.config`` instead."""
        warnings.warn(
            "ReadOFF.config is deprecated; use off.gmx.config instead.",
            DeprecationWarning, stacklevel=2
        )
        return self.gmx.config

    @config.setter
    def config(self, value):
        warnings.warn(
            "ReadOFF.config is deprecated; use off.gmx.config instead.",
            DeprecationWarning, stacklevel=2
        )
        self.gmx.config = value

    @property
    def nonbonded_tabpot(self):
        """Deprecated. Use ``off.gmx.nonbonded_tabpot`` instead."""
        warnings.warn(
            "ReadOFF.nonbonded_tabpot is deprecated; use off.gmx.nonbonded_tabpot instead.",
            DeprecationWarning, stacklevel=2
        )
        return self.gmx.nonbonded_tabpot

    @nonbonded_tabpot.setter
    def nonbonded_tabpot(self, value):
        warnings.warn(
            "ReadOFF.nonbonded_tabpot is deprecated; use off.gmx.nonbonded_tabpot instead.",
            DeprecationWarning, stacklevel=2
        )
        self.gmx.nonbonded_tabpot = value

    @property
    def bonded_tabpot(self):
        """Deprecated. Use ``off.gmx.bonded_tabpot`` instead."""
        warnings.warn(
            "ReadOFF.bonded_tabpot is deprecated; use off.gmx.bonded_tabpot instead.",
            DeprecationWarning, stacklevel=2
        )
        return self.gmx.bonded_tabpot

    @bonded_tabpot.setter
    def bonded_tabpot(self, value):
        warnings.warn(
            "ReadOFF.bonded_tabpot is deprecated; use off.gmx.bonded_tabpot instead.",
            DeprecationWarning, stacklevel=2
        )
        self.gmx.bonded_tabpot = value

    # ---------------------------------------------------------------------------
    # Deprecated wrappers — all GROMACS methods have moved to off.gmx
    # ---------------------------------------------------------------------------

    def gen_nonbonded_tabpot(self, *args, **kwargs):
        """Deprecated. Use ``off.gmx.gen_nonbonded_tabpot()`` instead."""
        warnings.warn(
            "ReadOFF.gen_nonbonded_tabpot() is deprecated; use off.gmx.gen_nonbonded_tabpot() instead.",
            DeprecationWarning, stacklevel=2
        )
        return self.gmx.gen_nonbonded_tabpot(*args, **kwargs)

    def gen_bonded_tabpot(self, *args, **kwargs):
        """Deprecated. Use ``off.gmx.gen_bonded_tabpot()`` instead."""
        warnings.warn(
            "ReadOFF.gen_bonded_tabpot() is deprecated; use off.gmx.gen_bonded_tabpot() instead.",
            DeprecationWarning, stacklevel=2
        )
        return self.gmx.gen_bonded_tabpot(*args, **kwargs)

    def gen_nonbonded_topology(self, *args, **kwargs):
        """Deprecated. Use ``off.gmx.gen_nonbonded_topology()`` instead."""
        warnings.warn(
            "ReadOFF.gen_nonbonded_topology() is deprecated; use off.gmx.gen_nonbonded_topology() instead.",
            DeprecationWarning, stacklevel=2
        )
        return self.gmx.gen_nonbonded_topology(*args, **kwargs)

    def gen_bonded_topology(self, *args, **kwargs):
        """Deprecated. Use ``off.gmx.gen_bonded_topology()`` instead."""
        warnings.warn(
            "ReadOFF.gen_bonded_topology() is deprecated; use off.gmx.gen_bonded_topology() instead.",
            DeprecationWarning, stacklevel=2
        )
        return self.gmx.gen_bonded_topology(*args, **kwargs)

    def write_nonbonded_tabpot(self, *args, **kwargs):
        """Deprecated. Use ``off.gmx.write_nonbonded_tabpot()`` instead."""
        warnings.warn(
            "ReadOFF.write_nonbonded_tabpot() is deprecated; use off.gmx.write_nonbonded_tabpot() instead.",
            DeprecationWarning, stacklevel=2
        )
        return self.gmx.write_nonbonded_tabpot(*args, **kwargs)

    def write_bonded_tabpot(self, *args, **kwargs):
        """Deprecated. Use ``off.gmx.write_bonded_tabpot()`` instead."""
        warnings.warn(
            "ReadOFF.write_bonded_tabpot() is deprecated; use off.gmx.write_bonded_tabpot() instead.",
            DeprecationWarning, stacklevel=2
        )
        return self.gmx.write_bonded_tabpot(*args, **kwargs)

    def set_config(self, **kwargs):
        """Deprecated. Use ``off.gmx.set_config()`` instead."""
        warnings.warn(
            "ReadOFF.set_config() is deprecated; use off.gmx.set_config() instead.",
            DeprecationWarning, stacklevel=2
        )
        return self.gmx.set_config(**kwargs)

    def get_config(self, key=None):
        """Deprecated. Use ``off.gmx.get_config()`` instead."""
        warnings.warn(
            "ReadOFF.get_config() is deprecated; use off.gmx.get_config() instead.",
            DeprecationWarning, stacklevel=2
        )
        return self.gmx.get_config(key)

    # ---------------------------------------------------------------------------
    # Non-GROMACS methods — remain on ReadOFF
    # ---------------------------------------------------------------------------

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
                atoms_to_mol = defaultdict(list)
                for key, value in self.charges.items():
                    atoms = [str(i) for i in value.keys()]
                    for atom in atoms:
                        atoms_to_mol[atom].append(key)
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
                        atomname, charge = parts[0], float(parts[1])
                        if current_mol is None:
                            print(f"Warning: Atom-charge pair '{line}' found before any molecule name. Adding atom to all possible molecules.")
                            mols_with_atom = atoms_to_mol[atomname]
                            for mol in mols_with_atom:
                                self.charges[mol][atomname] = charge

                        elif atomname not in self.charges[current_mol]:
                            print(f"Warning: Atom '{atomname}' not found in molecule '{current_mol}'. Skipping.")
                            continue

                        else:
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
