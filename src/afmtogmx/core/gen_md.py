import copy
from afmtogmx.core import functions, populate_bonded, residues
from afmtogmx.core.gmx_backend import GROMACSBackend
from afmtogmx.core.openmm_backend import OpenMMBackend
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
    openmm : OpenMMBackend
        OpenMM output backend. Provides ``gen_xml()`` for generating OpenMM
        force field XML files from the same parsed data.

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
        self.openmm = OpenMMBackend(self)

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

    def change_molecule(self, mol_name, reference_ff, atom_name_map=None, ref_mol_name=None):
        """Replace a molecule's force field parameters with a stored reference FF.

        Loads a reference `.off` file from the package's ``reference_ff/`` directory,
        replaces all bonded and water-water nonbonded parameters for ``mol_name`` with
        those from the reference, and renames all affected atom types using
        ``atom_name_map``. If a matching ``.charges`` file exists alongside the ``.off``
        file (same stem), it is loaded automatically.

        Cross-term pairs (one atom from ``mol_name``, one from another molecule) are
        preserved with their original fitted parameter values — only the atom type
        names on the ``mol_name`` side are renamed. This is correct because those
        cross terms were fitted against the reference FF as the fixed water model.

        Parameters
        ----------
        mol_name : str
            The molecule name as it appears in the ``.off`` file (e.g. ``'H2OQM'``).
        reference_ff : str
            Name of the reference FF file in ``src/afmtogmx/reference_ff/``, without
            the ``.off`` extension (e.g. ``'BLYPSP-4F'``).
        atom_name_map : dict, optional
            Maps atom type names in the current ``.off`` file to new names used in all
            output. Example: ``{'OW': 'OW_sp4f', 'HW': 'HW_sp4f'}``.
            Atom types not in the map are left unchanged. Defaults to ``{}``.
        ref_mol_name : str, optional
            Molecule name to use from the reference FF file. Required when the reference
            FF contains more than one molecule. If ``None`` and the reference FF has
            exactly one molecule, that molecule is used automatically.

        Returns
        -------
        ReadOFF
            Returns ``self`` for method chaining.

        Raises
        ------
        KeyError
            If ``mol_name`` is not found in ``self.bonded``, or ``ref_mol_name`` is
            specified but not found in the reference FF.
        FileNotFoundError
            If the reference FF file does not exist.
        ValueError
            If the reference FF file contains more than one molecule and
            ``ref_mol_name`` is not specified.

        Examples
        --------
        >>> off = ReadOFF('butanol_water.off')
        >>> off.change_molecule(
        ...     mol_name='H2OQM',
        ...     reference_ff='BLYPSP-4F',
        ...     ref_mol_name='H2OQM',
        ...     atom_name_map={'OW': 'OW_sp4f', 'HW': 'HW_sp4f'},
        ... )
        >>> # off.bonded['H2OQM'] now holds BLYPSP-4F bonded parameters
        >>> # off.nonbonded water-water pairs now use OW_sp4f / HW_sp4f atom types
        >>> # off.nonbonded cross pairs keep original values, water atoms renamed
        """
        from pathlib import Path

        if atom_name_map is None:
            atom_name_map = {}

        if mol_name not in self.bonded:
            raise KeyError(f"Molecule '{mol_name}' not found in self.bonded.")

        # Load the reference FF using the existing parser
        ref_path = Path(__file__).parent.parent / 'reference_ff' / f'{reference_ff}.off'
        ref = ReadOFF(str(ref_path))

        # Resolve which molecule to use from the reference FF
        if ref_mol_name is None:
            if len(ref.bonded) == 1:
                ref_mol_name = list(ref.bonded.keys())[0]
            else:
                raise ValueError(
                    f"Reference FF '{reference_ff}' contains multiple molecules "
                    f"({list(ref.bonded.keys())}). Specify ref_mol_name."
                )
        elif ref_mol_name not in ref.bonded:
            raise KeyError(
                f"Molecule '{ref_mol_name}' not found in reference FF '{reference_ff}'."
            )

        # Auto-load a matching .charges file if one exists alongside the .off file
        ref_charges_path = ref_path.with_suffix('.charges')
        if ref_charges_path.exists():
            ref.load_charges_from_file(str(ref_charges_path))

        # Atom types that belong to mol_name (excluding NETF/TORQ)
        mol_atoms = {
            entry[1]
            for entry in self.bonded[mol_name]['ATO']['All'].values()
            if entry[1] not in ('NETF', 'TORQ')
        }

        # Atom types that belong to ref_mol_name in the reference FF
        ref_mol_atoms = {
            entry[1]
            for entry in ref.bonded[ref_mol_name]['ATO']['All'].values()
            if entry[1] not in ('NETF', 'TORQ')
        }

        # Replace bonded parameters
        self.bonded[mol_name] = ref.bonded[ref_mol_name]

        # Replace charges — rename keys via atom_name_map; look up values by old name
        # (before renaming) since ref_charges is indexed by the reference FF's atom type names
        ref_charges = ref.charges.get(ref_mol_name, {})
        self.charges[mol_name] = {
            atom_name_map.get(old, old): ref_charges.get(old, 0.0)
            for old in self.charges.get(mol_name, {})
        }

        # Update nonbonded pairs:
        #   - pure mol_name pairs  → remove (reference FF supplies replacements)
        #   - cross pairs          → rename mol_name atom type, keep parameter values
        pairs_to_remove = []
        cross_pairs = {}
        for pair, params in self.nonbonded.items():
            a1, a2 = pair
            a1_in = a1 in mol_atoms
            a2_in = a2 in mol_atoms
            if a1_in and a2_in:
                pairs_to_remove.append(pair)
            elif a1_in:
                pairs_to_remove.append(pair)
                cross_pairs[(atom_name_map.get(a1, a1), a2)] = params
            elif a2_in:
                pairs_to_remove.append(pair)
                cross_pairs[(a1, atom_name_map.get(a2, a2))] = params

        for pair in pairs_to_remove:
            del self.nonbonded[pair]
        self.nonbonded.update(cross_pairs)

        # Add reference FF nonbonded pairs for ref_mol_name only, applying atom_name_map
        renamed_ref_nonbonded = {
            (atom_name_map.get(a1, a1), atom_name_map.get(a2, a2)): params
            for (a1, a2), params in ref.nonbonded.items()
            if a1 in ref_mol_atoms and a2 in ref_mol_atoms
        }
        self.nonbonded.update(renamed_ref_nonbonded)

        # Rebuild residues from updated bonded data
        self.residues = {
            "Definitions": {
                k: {'All': functions._remove_netf_torq_atname(v['ATO']['All'])}
                for k, v in self.bonded.items()
            },
            "Residues": {
                k: {'All': [functions._remove_netf_torq_atnum(v['ATO']['All'])]}
                for k, v in self.bonded.items()
            },
        }

        return self

    def populate_bonded(self, directory, new_mol_name, remove_molecules=None):
        """Build a new ReadOFF whose bonded section for ``new_mol_name`` is
        reconstructed from a user-supplied topology directory.

        The directory must contain ``atoms.dat``, ``bonds.dat``,
        ``valid_dihedrals.dat``, and ``parameters.dat``. See
        :mod:`afmtogmx.core.populate_bonded` for the file formats.

        Co-fitted molecules in the original ``self`` (e.g. a co-fitted water
        model) are carried through unchanged. ``remove_molecules`` lets the
        caller drop the original monomer being replaced. ``self.nonbonded``
        is filtered to atom-type pairs whose types still appear in some
        surviving molecule.

        Parameters
        ----------
        directory : str
            Path to the topology-input directory.
        new_mol_name : str
            Name for the assembled molecule in the returned ReadOFF.
        remove_molecules : list of str, optional
            Molecule names in the current ReadOFF to drop from the new object.

        Returns
        -------
        ReadOFF
            A new ReadOFF instance. ``self`` is not modified.

        Raises
        ------
        ValueError
            If a required file is missing, an atom type in ``atoms.dat`` is
            not present in any surviving molecule, or any bond/angle/dihedral
            signature cannot be resolved against ``parameters.dat``.
        """
        if remove_molecules is None:
            remove_molecules = []

        new = copy.deepcopy(self)

        # Atom types known to the original FF: union of ATO types across every
        # molecule (including those about to be removed, since their nonbonded
        # entries carry across as long as some surviving molecule uses the same
        # type — typical when CB7 reuses C7F's types).
        atom_type_universe = set()
        for mol_data in new.bonded.values():
            for entry in mol_data['ATO']['All'].values():
                atype = entry[1]
                if atype not in ('NETF', 'TORQ'):
                    atom_type_universe.add(atype)

        for name in remove_molecules:
            new.bonded.pop(name, None)
            new.charges.pop(name, None)
            new.residues['Definitions'].pop(name, None)
            new.residues['Residues'].pop(name, None)

        new.bonded[new_mol_name] = populate_bonded.build_new_molecule_bonded(
            directory, atom_type_universe
        )

        # Charges: one entry per atom type in the new molecule, default 0.0
        # (matching the __init__ convention of keying by atom type via pair[1]).
        new_types = {entry[1] for entry in new.bonded[new_mol_name]['ATO']['All'].values()}
        new.charges[new_mol_name] = {atype: 0.0 for atype in new_types}

        # Prune nonbonded to atom-type pairs whose types still appear somewhere.
        surviving_types = set()
        for mol_data in new.bonded.values():
            for entry in mol_data['ATO']['All'].values():
                atype = entry[1]
                if atype not in ('NETF', 'TORQ'):
                    surviving_types.add(atype)
        new.nonbonded = {
            pair: params
            for pair, params in new.nonbonded.items()
            if pair[0] in surviving_types and pair[1] in surviving_types
        }

        # Add residues entry for new_mol_name, matching __init__'s one-liner.
        ato_all = new.bonded[new_mol_name]['ATO']['All']
        new.residues['Definitions'][new_mol_name] = {
            'All': functions._remove_netf_torq_atname(ato_all)
        }
        new.residues['Residues'][new_mol_name] = {
            'All': [functions._remove_netf_torq_atnum(ato_all)]
        }

        return new

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
