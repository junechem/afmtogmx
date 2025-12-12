"""This module contains functions mainly related to the production of topology files; for now 
this only generates .top files
"""
import re


def _find_keyword_location(template_file, keyword="", begin=0):
    """Find the beginning and end character locations of a keyword section in a template file.

    This function searches for a specific keyword section (e.g., `[ nonbond_params ]`)
    within a GROMACS topology template file string. It returns the character indices
    of the keyword's start and the subsequent blank line, while accounting for comments.

    Parameters
    ----------
    template_file : str
        The full content of the topology template file as a string.
    keyword : str, optional
        The keyword to search for (e.g., "nonbond_params"). Defaults to an empty string.
    begin : int, optional
        The starting character index in `template_file` to begin the search. Defaults to 0.

    Returns
    -------
    tuple of (int, int)
        A tuple containing the starting and ending character indices of the found
        keyword section.

    Raises
    ------
    SystemExit
        If the specified `keyword` cannot be found in the `template_file` after
        the `begin` index, the program exits with an error message.
    """

def _gen_nonbonded_string(scale_C6, special_pairs, name_translation, nonbonded):
    """Generate the `[ nonbond_params ]` section as a string.

    This function iterates through the filtered nonbonded interactions to
    determine the appropriate C6 and C12 parameters for each atom pair.
    It then formats these parameters into a string suitable for inclusion
    in the `[ nonbond_params ]` section of a GROMACS topology file.

    Parameters
    ----------
    scale_C6 : bool
        If `True`, C6 parameters will be adjusted to allow for correct
        dispersion corrections in GROMACS simulations.
    special_pairs : dict
        A dictionary defining custom attractive interactions for specific
        pairs.
    name_translation : dict
        A dictionary for translating atom names.
    nonbonded : dict
        A filtered dictionary of nonbonded interactions, typically from
        `_filter_nonbonded`.

    Returns
    -------
    str
        A formatted string representing the `[ nonbond_params ]` section.

    Raises
    ------
    SystemExit
        If `special_pairs` is used in conjunction with `scale_C6 = True`,
        as this combination is not supported due to potential conflicts
        in parameter scaling.
        If more than one attractive interaction is found for a pair and
        `scale_C6 = True`, as this also leads to conflicts.
    """

    if special_pairs and scale_C6:
        print(
            "Using special_pairs when also scaling C6 interactions is not supported. Please re-run with scale_C6 = False")
        exit(1)
    default_attractive = ['POW_6', 'DPO_6', 'SRD_6']

    nonbonded_string = ""

    for pair, int_dict in nonbonded.items():  # determine what the C6, C12 should be for each At1~At2 nonbonded interaction
        num_attr = 0
        if pair not in special_pairs:
            C6, C12 = 0.0, 0.0
            for interaction in int_dict:
                if interaction == 'BUC':  # Special case for Buckingham potential: U(r) = A*exp(-B*r) - C/r^6
                    for param_set in int_dict[interaction]:
                        num_attr += 1
                        if num_attr > 1 and scale_C6:
                            print(
                                "Cannot have greater than 1 attractive interactions and properly scale C6 for dispersion "
                                "corrections. Set 'scale_C6' to false and restart to continue generating nonbonded")
                            exit(1)
                        # BUCK third parameter is the C6 coefficient (attractive r^-6 term)
                        C6 = abs(param_set[1] * 4.184 * 1E-6)
                        C12 = 1  # Also has repulsive exponential component
                elif interaction not in default_attractive and interaction != 'THC':
                    C12 = 1
                elif interaction == 'THC':
                    print("UNABLE TO HANDLE THC FUNCTIONS AT THIS TIME: EDIT CODE YOURSELF OR ASK RAY TO DO IT")
                else:
                    for param_set in int_dict[interaction]:
                        num_attr += 1
                        if num_attr > 1 and scale_C6:
                            print(
                                "Cannot have greater than 1 attractive interactions and properly scale C6 for dispersion "
                                "corrections. Set 'scale_C6' to false and restart to continue generating nonbonded")
                            exit(1)
                        C6 = abs(param_set[0] * 4.184 * 1E-6)

        else:
            C6, C12 = 1.0, 1.0

        nonbonded_string += single_nonbonded_pair_string(pair, name_translation, C6, C12)

    return nonbonded_string


def single_nonbonded_pair_string(pair, name_translation, C6, C12):
    """Generate a formatted string for a single nonbonded atom pair.

    This function creates a line of text formatted according to GROMACS
    `[ nonbond_params ]` section requirements for a given atom pair,
    including their C6 and C12 parameters. Atom names can be translated.

    Parameters
    ----------
    pair : tuple of str
        A tuple containing the atom names for the nonbonded pair, e.g., `('C', 'H')`.
    name_translation : dict
        A dictionary for translating atom names. Format:
        `{'.offAtom1': '.topAtom1', ...}`. If an atom is not found in
        `name_translation`, its original name is used.
    C6 : float
        The C6 parameter for the Lennard-Jones potential.
    C12 : float
        The C12 parameter for the Lennard-Jones potential.

    Returns
    -------
    str
        A string formatted for the `[ nonbond_params ]` section,
        e.g., `'At1    At2    1      C6_value           C12_value\\n'`.
    """
    if pair[0] in name_translation:
        At1 = name_translation[pair[0]]
    else:
        At1 = pair[0]
    if pair[1] in name_translation:
        At2 = name_translation[pair[1]]
    else:
        At2 = pair[1]

    return f'{At1:<7}{At2:<7}{1:<7}{C6:<20.12E}{C12:<20.12E}\n'


def _gen_bonded_section_strings(name_translation, filtered_bonded, charges, molname, bonded_tabpot):
    """Generate a dictionary of formatted strings for all bonded sections of a molecule.

    This function processes the filtered bonded interaction data for a given
    molecule and generates formatted string content for each relevant GROMACS
    topology keyword section (e.g., `[ atoms ]`, `[ bonds ]`, `[ angles ]`).
    It populates a skeleton dictionary with these strings.

    Parameters
    ----------
    name_translation : dict
        A dictionary for translating atom names, e.g., `{'.offAtom1' : '.topAtom1', ...}`.
    filtered_bonded : dict
        A filtered dictionary of bonded interactions for a molecule, with
        empty interactions removed, typically from `functions._filter_bonded`.
    charges : dict
        The `self.charges` dictionary from the `ReadOFF` class, containing
        atomic charges.
    molname : str
        The name of the molecule being processed.
    bonded_tabpot : dict
        A dictionary of bonded tabulated potentials, required for `QUA` and
        `QBB` interactions.

    Returns
    -------
    dict
        A dictionary where keys are GROMACS topology keywords (e.g., 'atoms',
        'bonds', 'angles') and values are the formatted string content for
        those sections for the specified molecule.
    """
    bonded_skeleton = _gen_bonded_skeleton(filtered_bonded)

    for interaction, int_dict in filtered_bonded.items():
        if interaction == 'ATO':
            bonded_skeleton['atoms'] += _gen_bonded_atoms(name_translation, charges, int_dict, molname)
            ### Below should be added to its own function, I am just lazy right now
            if 'Virtual' in int_dict:
                if 'virtual_sitesn' in bonded_skeleton:
                    pass
                else:
                    bonded_skeleton['virtual_sitesn'] = ""
                for atoms, params in int_dict['Virtual'].items():
                    param_list = list(params[1:])
                    param_list = [i for i in param_list if i != "+"]
                    param_list = param_list[::-1]
                    param_string = f'{atoms[0]:8}{3:8}'
                    for item in param_list:
                        param_string += f'{item:>5}'
                        param_string += "  "
                    param_string = re.sub('\(', ' ', param_string)
                    param_string = re.sub('\)', ' ', param_string)
                    param_string += "\n"
                    bonded_skeleton['virtual_sitesn'] += param_string
        if interaction == 'BON':
            bonded_skeleton['bonds'] += _gen_bonded_bonds(int_dict, molname, bonded_tabpot)
        if interaction == 'ANG':
            bonded_skeleton['angles'] += _gen_bonded_angles(int_dict)
        if interaction == 'BD3':
            bonded_skeleton['bonds'] += _gen_bonded_bd3(int_dict, molname=molname, bonded_tabpot=bonded_tabpot,
                                                        bonds=True)
            bonded_skeleton['angles'] += _gen_bonded_bd3(int_dict, bonds=False)
        if interaction == 'DIH':
            bonded_skeleton['dihedrals'] += _gen_dihedrals(int_dict)
        if interaction == 'CDI':
            bonded_skeleton['cmap'] += _gen_cmap(int_dict)
        if interaction == 'EXC':
            bonded_skeleton['exclusions'] += _gen_exclusions(int_dict)
    return bonded_skeleton


def _gen_bonded_skeleton(filtered_bonded):
    """Generate an empty skeleton dictionary for bonded topology sections.

    This function initializes a dictionary with appropriate GROMACS topology
    keyword keys (e.g., 'atoms', 'bonds', 'angles') based on the types of
    bonded interactions present in `filtered_bonded`. Each key is mapped
    to an empty string, which will later be populated with formatted data.

    Parameters
    ----------
    filtered_bonded : dict
        A dictionary of bonded interactions for a single molecule, with
        empty interactions removed, typically from `functions._filter_bonded`.

    Returns
    -------
    dict
        An initialized dictionary where keys are GROMACS topology keywords
        and values are empty strings.
        Example: `{'atoms': '', 'bonds': '', 'angles': '', ...}`
    """
    bonded_skeleton = dict()
    translated_names = {'ATO': ['atoms'], 'BON': ['bonds'], 'ANG': ['angles'], 'BD3': ['bonds', 'angles'],
                        'DIH': ['dihedrals'], 'CDI': ['cmap'], 'EXC': ['exclusions']}
    for interaction in filtered_bonded:
        for keyword in translated_names[interaction]:
            bonded_skeleton[keyword] = ""  # generate key for each match in translated_names
    return bonded_skeleton


def _gen_bonded_atoms(name_translation, charges, int_dict, molname):
    """Generate the formatted string for the `[ atoms ]` section of a GROMACS topology.

    This function takes atom definitions, charges, and name translation data
    to produce a formatted string representing the `[ atoms ]` section for
    a specific molecule in a GROMACS topology file.

    Parameters
    ----------
    name_translation : dict
        A dictionary for translating atom names.
    charges : dict
        The `self.charges` dictionary containing atomic charges for each molecule.
    int_dict : dict
        The interaction dictionary for the 'ATO' (atom) section of a molecule.
    molname : str
        The name of the molecule.

    Returns
    -------
    str
        A formatted string content for the `[ atoms ]` section, including
        atom number, atom name, residue number, residue name, atom type,
        charge group, and charge.
    """
    resnum = 1  # To Do: allow for dictionary to delineate which atoms go to which residue, have more than just 1 residue
    resname = molname
    bonded_atoms_string = ""
    for atnum, atname_pair in int_dict['All'].items():
        if atname_pair[0] in name_translation:
            atom_name = name_translation[atname_pair[0]]
            charge = charges[molname][atname_pair[1]]
        elif 'NETF' not in atname_pair and 'TORQ' not in atname_pair:
            atom_name = atname_pair[0]
            charge = charges[molname][atname_pair[1]]
        else:
            continue
        bonded_atoms_string += f"{int(atnum):<8}{atom_name:<8}{resnum:<8}{resname:<8}{atom_name:<8}{int(atnum):<8}{charge:<8.5f}\n"
    return bonded_atoms_string


def _gen_bonded_bonds(int_dict, molname, bonded_tabpot):
    """Generate the formatted string for the `[ bonds ]` section of a GROMACS topology.

    This function processes bond interaction data for a given molecule and
    generates a formatted string suitable for the `[ bonds ]` section. It
    handles both harmonic (`HAR`) and quartic tabulated (`QUA`) bond types,
    referencing `bonded_tabpot` for `QUA` bonds.

    Parameters
    ----------
    int_dict : dict
        The interaction dictionary for the 'BON' (bonds) section of a molecule.
    molname : str
        The name of the molecule.
    bonded_tabpot : dict
        A dictionary of bonded tabulated potentials, required for `QUA` bonds.

    Returns
    -------
    str
        A formatted string content for the `[ bonds ]` section.

    Raises
    ------
    SystemExit
        If a `QUA` bond is defined but its corresponding `molname` or parameters
        are not found in `bonded_tabpot`, indicating a potential mismatch
        between `gen_bonded_tabpot()` and `gen_bonded_topology()`.
    """
    bonded_bonds_string = ""

    for bond_type, bond_type_dict in int_dict.items():
        if bond_type == 'HAR':
            for params, atnum_list in bond_type_dict.items():
                P1, P2 = 1E-1 * params[0], 4.184E2 * params[1]
                for at1, at2 in atnum_list:
                    bonded_bonds_string += f"{at1:8}{at2:8}{1:8}{P1:10.5f}{P2:20.2f}\n"
        if bond_type == 'QUA':
            for params, atnum_list in bond_type_dict.items():
                try:
                    table_number = bonded_tabpot[molname][params][0]
                except KeyError as e:
                    print('ERROR with gen_bonded_topology()')
                    print(f'There is no {molname} found in your bonded_tabpot file. This likely means that you used'
                          f' incl_mol when calling gen_bonded_tabpot() and not when using gen_bonded_topology. Be sure'
                          f' to pass incl_mol to gen_bonded_topology() if it was used for previous functions.')
                    exit(1)
                for at1, at2 in atnum_list:
                    bonded_bonds_string += f"{at1:8}{at2:8}{8:8}{int(table_number):8}{1.000:8}\n"
    return bonded_bonds_string


def _gen_bonded_angles(int_dict):
    """Generate the formatted string for the `[ angles ]` section of a GROMACS topology.

    This function processes angle interaction data for a given molecule and
    generates a formatted string suitable for the `[ angles ]` section. It
    handles both harmonic (`HAR`) and quartic (`QUA`) angle types.

    Parameters
    ----------
    int_dict : dict
        The interaction dictionary for the 'ANG' (angles) section of a molecule.

    Returns
    -------
    str
        A formatted string content for the `[ angles ]` section.
    """
    bonded_angles_string = ""

    for angle_type, angle_type_dict in int_dict.items():
        if angle_type == 'HAR':
            for params, atnum_list in angle_type_dict.items():
                P1, P2 = params[0], 4.184 * params[1]
                for at1, at2, at3 in atnum_list:
                    bonded_angles_string += f"{at1:8}{at2:8}{at3:8}{1:8}{P1:10.5f}{P2:10.4f}\n"
        elif angle_type == "QUA":
            for params, atnum_list in angle_type_dict.items():
                P1, P2, P3, P4 = params[0], 4.184 * params[1], 4.184 * params[2], 4.184 * params[3]
                for at1, at2, at3 in atnum_list:
                    bonded_angles_string += f"{at1:8}{at2:8}{at3:8}{6:8}{P1:10.5f}{0.0:8}{P2:10.4f}{P3:10.4f}{P4:10.4f}"

    return bonded_angles_string


def _gen_bonded_bd3(int_dict, molname="", bonded_tabpot={}, bonds=False):
    """Generate formatted strings for BD3 (bond-dipole 3-center) interactions.

    This function processes BD3 interaction data for a given molecule and
    generates a formatted string suitable for either the `[ bonds ]` or
    `[ angles ]` section, depending on the `bonds` parameter. It handles
    'QBB' (quartic bond-bond) and 'MUB' (Morse Urey-Bradley) interaction types.

    Parameters
    ----------
    int_dict : dict
        The interaction dictionary for the 'BD3' section of a molecule.
    molname : str, optional
        The name of the molecule. Required for 'QBB' interactions when looking
        up `bonded_tabpot`. Defaults to an empty string.
    bonded_tabpot : dict, optional
        A dictionary of bonded tabulated potentials, used for 'QBB' interactions.
        Defaults to an empty dictionary.
    bonds : bool, optional
        If `True`, the function generates a string for the `[ bonds ]` section
        (specifically for 'QBB' interactions). If `False`, it generates a
        string for the `[ angles ]` section. Defaults to `False`.

    Returns
    -------
    str
        A formatted string content for either the `[ bonds ]` or `[ angles ]`
        section, representing the BD3 interactions.
    """
    bonded_bd3_string = ""
    for bd3_type, bd3_type_dict in int_dict.items():
        if bd3_type == 'QBB' and bonds:
            for params, atnum_list in bd3_type_dict.items():
                table_number = bonded_tabpot[molname][params][0]
                for at1, at2, at3 in atnum_list:
                    a1, a2 = sorted([at1, at2])
                    a3, a4 = sorted([at2, at3])
                    bonded_bd3_string += f"{a1:8}{a2:8}{8:8}{int(table_number):8}{1.000:8}\n"
                    bonded_bd3_string += f"{a3:8}{a4:8}{8:8}{int(table_number):8}{1.000:8}\n"
        elif bd3_type == 'QBB' and not bonds:
            for params, atnum_list in bd3_type_dict.items():
                for at1, at2, at3 in atnum_list:
                    P1, P2, P3 = 1E-1 * params[0], 1E-1 * params[0], 4.184E2 * params[1]
                    bonded_bd3_string += f"{at1:8}{at2:8}{at3:8}{3:8}{P1:10.5f}{P2:10.5f}{P3:20.2f}\n"
        elif bd3_type == "MUB":
            for params, atnum_list in bd3_type_dict.items():
                for at1, at2, at3 in atnum_list:
                    P1, P2, P3, P4 = 4.184E2 * params[0], 1E-1 * params[1], 1E-1 * params[2], 1E-1 * params[3]
                    print("DUE TO MANUAL BEING UNCLEAR< NOT SURE IF BONDED 'MUB' INTERACTION HAS PROPER PARAMETERS")
                    bonded_bd3_string += f"{at1:8}{at2:8}{at3:8}{4:8}{P2:10.5f}{P3:10.5F}{P4:10.5f}{P1:10.5f}\n"

    return bonded_bd3_string


def _gen_dihedrals(int_dict):
    """Generate the formatted string for the `[ dihedrals ]` section of a GROMACS topology.

    This function processes dihedral interaction data for a given molecule and
    generates a formatted string suitable for the `[ dihedrals ]` section. It
    handles harmonic (`HAR`), Ryckaert-Bellemans (`NCO`), and cosine-series (`COS`)
    dihedral types.

    Parameters
    ----------
    int_dict : dict
        The interaction dictionary for the 'DIH' (dihedrals) section of a molecule.

    Returns
    -------
    str
        A formatted string content for the `[ dihedrals ]` section.
    """
    bonded_dihedrals_string = ""
    for dihedral_type, dihedral_type_dict in int_dict.items():
        if dihedral_type == 'HAR':
            for params, atnum_list in dihedral_type_dict.items():
                P1, P2 = params[0], 4.184 * params[1]
                for at1, at2, at3, at4 in atnum_list:
                    bonded_dihedrals_string += f"{at1:8}{at2:8}{at3:8}{at4:8}{2:8}{P2:10.4f}{P1:10.4f}\n"
        if dihedral_type == 'NCO':
            for params, atnum_list in dihedral_type_dict.items():
                P1, P2, P3 = 4.184 * params[0], params[1], params[2]
                for at1, at2, at3, at4 in atnum_list:
                    bonded_dihedrals_string += f"{at1:8}{at2:8}{at3:8}{at4:8}{9:8}{P3:10.5f}{P1:10.5f}{P2:8}\n"
        if dihedral_type == 'COS':
            for params, atnum_list in dihedral_type_dict.items():
                P1, P2, P3 = 4.184 * params[0], params[1], params[2]
                for at1, at2, at3, at4 in atnum_list:
                    bonded_dihedrals_string += f"{at1:8}{at2:8}{at3:8}{at4:8}{9:8}{P3:10.5f}{P1:10.5f}{P2:8}\n"

    return bonded_dihedrals_string


def _gen_cmap(int_dict):
    """Generate the formatted string for the `[ cmap ]` section of a GROMACS topology.

    This function processes cmap (cross-molecular angle potential) interaction
    data for a given molecule and generates a formatted string suitable for
    the `[ cmap ]` section.

    Parameters
    ----------
    int_dict : dict
        The interaction dictionary for the 'CDI' (cmap dihedrals) section
        of a molecule.

    Returns
    -------
    str
        A formatted string content for the `[ cmap ]` section.
    """
    cmap_string = ""
    for cmap_type, cmap_type_dict in int_dict.items():
        for atnum_list in cmap_type_dict.items():
            for at1, at2, at3, at4, at5 in atnum_list:
                cmap_string += f"{at1:8}{at2:8}{at3:8}{at4:8}{at5:8}\n"
    return cmap_string


def _gen_exclusions(int_dict):
    """Generate the formatted string for the `[ exclusions ]` section of a GROMACS topology.

    This function processes exclusion interaction data for a given molecule and
    generates a formatted string suitable for the `[ exclusions ]` section.

    Parameters
    ----------
    int_dict : dict
        The interaction dictionary for the 'EXC' (exclusions) section of a molecule.

    Returns
    -------
    str
        A formatted string content for the `[ exclusions ]` section.
    """
    exclusions_string = ""
    for atom_list in int_dict:
        for atom in atom_list:
            if int(atom) != 0:
                exclusions_string += f'{atom:5}'
        exclusions_string += '\n'
    return exclusions_string


def _find_moleculetypes(topology):
    """Find all `[ moleculetype ]` sections within a topology string, accounting for comments.

    This function searches for `[ moleculetype ]` declarations in a GROMACS
    topology file's content, identifying their start and end character locations.
    It intelligently ignores sections that are commented out.

    Parameters
    ----------
    topology : str
        The full content of the topology template file as a string.

    Returns
    -------
    list of tuples
        A list of tuples, where each tuple contains:
        - `(begin_match, end_match)`: A tuple representing the start and end
          character indices of the `[ moleculetype ]` section in the `topology` string.
        - `matching_string`: The uncommented content of the matched section.
    """
    uncommented_locations = []
    matches = re.finditer(r'\[[\t ]*moleculetype[\t ]*][\s\S]*?^\s*?$', topology, re.MULTILINE)
    for match in matches:
        beg_match = match.span()[0]
        prev_newline = topology[:beg_match].rfind('\n')

        if ';' in topology[prev_newline:beg_match]:
            pass
        else:
            uncommented_match = re.sub(';.*', '', topology[match.span()[0]:match.span()[1]])
            match = (match.span(), uncommented_match)
            uncommented_locations.append(match)
    return uncommented_locations


def _gen_bonded_string(topology, locations, topology_strings_dict):
    """Generate the complete, formatted GROMACS topology file as a string.

    This function reconstructs the entire topology file by combining the
    header from the template, the processed `[ moleculetype ]` sections,
    and any other sections that need to be preserved.

    Parameters
    ----------
    topology : str
        The original content of the topology template file as a string.
    locations : list
        A list of `[ moleculetype ]` section locations, typically the output
        of `_find_moleculetypes()`.
    topology_strings_dict : dict
        A dictionary containing the generated strings for each molecule's
        bonded sections, typically the output of `_gen_bonded_section_strings()`.
        Format: `{molname: {'.top_keyword' : section_string}, ...}`.

    Returns
    -------
    str
        The complete and properly formatted GROMACS topology file content as a string.
    """
    header_string = topology[:locations[0][0][0]]
    body_string = ""
    counter = 0
    for moleculetype in locations:
        body_string += topology[moleculetype[0][0]:moleculetype[0][1]]
        # Below code before the 'for' loop ensures that items below [moleculetypes] are preserved in the final topology
        # file
        next_keyword = re.search(r'\[[\t ]*(\w*)[\t ]*][\s\S]*?^\s*?$', topology[moleculetype[0][1]:], re.MULTILINE)
        beg_next_keyword = next_keyword.span()[0]
        prev_newline = moleculetype[0][1] + topology[moleculetype[0][1]: moleculetype[0][1] + beg_next_keyword].rfind('\n')
        begin, end = moleculetype[0][1], prev_newline
        body_string += topology[begin:end]
        for molname in topology_strings_dict:
            in_off = False
            if molname in moleculetype[1]:
                in_off = True
                body_string += _gen_molname_bonded(topology=topology[moleculetype[0][1]:],
                                                   bonded=topology_strings_dict[molname])
                break
        if in_off:  # write out rest of information below [ moleculetype ] section
            pass
        else:
            if moleculetype == locations[-1]:
                body_string += topology[moleculetype[0][1]:]
            else:
                body_string += topology[moleculetype[0][1]:locations[counter + 1][0][0]]
        counter += 1

    total_topology = header_string + body_string

    return total_topology


def _gen_molname_bonded(topology, bonded):
    """Generate the properly formatted `[ moleculetype ]` section for a single molecule.

    This function takes the template topology content and the pre-generated
    bonded section strings for a specific molecule, then inserts these strings
    into the correct locations within the `[ moleculetype ]` section.

    Parameters
    ----------
    topology : str
        A segment of the template topology file string, starting from the
        `[ moleculetype ]` keyword.
    bonded : dict
        A dictionary containing the generated strings for a molecule's bonded
        sections. Format: `{'.top_keyword' : section_string}, ...`.

    Returns
    -------
    str
        The complete and properly formatted `[ moleculetype ]` section
        content for the specified molecule as a string.
    """
    body_str = ""
    matches = [item for item in re.finditer(r'\[[\t ]*(\w*)[\t ]*][\s\S]*?^\s*?$', topology, re.MULTILINE)]
    counter = 0
    for match in matches:
        keyword = match.group(1)
        last_match = False
        if match == matches[-1]:
            last_match = True
        beg_match = match.span()[0]
        prev_newline = topology[:beg_match].rfind('\n')
        if ';' in topology[prev_newline:beg_match] and last_match:
            body_str += topology[beg_match:]
            continue
        elif ';' in topology[prev_newline:beg_match]:
            body_str += topology[prev_newline:matches[counter + 1].span()[0]]
        elif keyword in bonded:
            body_str += topology[match.span()[0]: match.span()[1]]
            body_str += bonded[keyword]
            if not last_match:
                body_str += topology[match.span()[1]:matches[counter + 1].span()[0]]
            else:
                body_str += topology[match.span()[1]:]
        elif keyword == 'system' or keyword == 'molecules':
            body_str += topology[match.span()[0]:]
            break
        elif keyword == 'moleculetype':
            break
        else:
            if last_match:
                body_str += topology[match.span()[0]:]
            else:
                body_str += topology[match.span()[0]:matches[counter + 1].span()[0]]

        counter += 1

    return body_str


def _write_topology(final_topology, write_to):
    """Write the final GROMACS topology string to a specified file.

    Parameters
    ----------
    final_topology : str
        The complete, formatted GROMACS topology file content as a string.
    write_to : str
        The path to the output `.top` file, including the filename.
    """
    with open(f'{write_to}', 'w') as f:
        f.write(final_topology)
