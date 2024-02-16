import re


def _find_keyword_location(template_file, keyword="", begin=0):
    """Returns beginning (location of initial \[) and ending (location of next blank line)"""

    # This script could be improved upon/combined with scripts for matching keywords further down
    found = False
    expression = rf'\[[\t ]*{keyword}[\t ]*][\s\S]*?^\s*?$'  # painful regex to match up to next blank line
    matches = re.finditer(expression, template_file, re.MULTILINE)  # find all matches
    for match in matches:  # check to get first match that has location that is greater than begin
        if match.span()[0] >= begin:
            poss_beg = match.span()[0]
            prev_newline = template_file[:poss_beg].rfind('\n')
            if prev_newline == poss_beg:  # easy check if '[' is at beginning of line
                found = True
                break
            else:
                if ";" in template_file[prev_newline:poss_beg]:  # check for comments
                    pass
                else:
                    found = True  # if no comments, we found the match
                    break
        else:
            pass  # skip characters with location < begin

    if not found:  # If a match could not be found at all
        print(
            f"Could not find valid KEYWORD '{keyword}' between character at LOCATION '{begin}' and the end of the file."
            "Please check your inputs")
        exit(1)
    else:  # If a match was found, return beginning and ending of match
        return match.span()[0], match.span()[1]


def _gen_nonbonded_string(scale_C6, special_pairs, name_translation, nonbonded):
    """Returns completed [ nonbonded_params ] section as a string. Atoms in special_pairs will have C6, C12 set to 1.

    :param scale_C6: Boolean, checking whether C6 scaling for dispersion corrections is on. Incompatible with special_pairs
    :param special_pairs: Dictionary, same as gen_nonbonded_tabpot
    :param name_translation: Dictionary, containing {'AtInOff' : 'AtIn.top', ...}. If atom not found in name_translation, name from .off is used
    :param nonbonded: filtered self.nonbonded
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
                if interaction not in default_attractive and interaction != 'THC':
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
    """Returns string for given pair, C6, C12, and name_translation, that will be written under [ nonbonded_params ] section

    :param pair: tuple of atom pair
    :param name_translation: dictionary with {'.offAtom1': '.topAtom1',...}
    :param C6: float, C6 parameter
    :param C12: float, C12 parameter
    :return: string with proper formatting for nonbonded_params section
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
    """Generates dictionary containing topology keywords Key1, Key2,..., with values Key1Section, Key2Section,..., as key : value pairs in a dictionary

    :param name_translation: dictionary, {'.offAtom1' : '.topAtom1',...}
    :param filtered_bonded: self.bonded with empty interactions removed
    :param charges: self.charges dictionary
    :param molname: molname string
    :param bonded_tabpot: dictionary, bonded_tabulated parameters. Generally empty
    :return: dictionary, containing necessary keywords and strings that should be written in those sections for a molname in a topology file
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
    """Takes in filtered_bonded (bonded dict with empty items removed) and outputs skeleton with proper format, excluding virtual sites

    :param filtered_bonded: dictionary, output of functions._filter_bonded
    :return: dictionary with format {'keyword' : ''}
    """
    bonded_skeleton = dict()
    translated_names = {'ATO': ['atoms'], 'BON': ['bonds'], 'ANG': ['angles'], 'BD3': ['bonds', 'angles'],
                        'DIH': ['dihedrals'], 'CDI': ['cmap'], 'EXC': ['exclusions']}
    for interaction in filtered_bonded:
        for keyword in translated_names[interaction]:
            bonded_skeleton[keyword] = ""  # generate key for each match in translated_names
    return bonded_skeleton


def _gen_bonded_atoms(name_translation, charges, int_dict, molname):
    """Generate 'atoms' section string for moleculetype 'molname'

    :param name_translation: dictionary
    :param charges: self.charges
    :param int_dict: dictionary, for interaction 'ATO'
    :param molname: molname
    :return: string for 'atoms' section in topology file of moleculetype 'molname'
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
    """Generate 'bonds' section string for moleculetype 'molname'

    :param int_dict: dictionary, 'BON' interaction dictionary
    :param molname: str, molname
    :param bonded_tabpot: dictionary, bonded tabulated potentials (if necessary)
    :return: string for 'bonds' section in topology file of moleculetype 'molname'
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
    """Generate 'angles' section string for moleculetype 'molname'

    :param int_dict: dictionary,'ANG' interaction dictionary for molname
    :return: string for 'angles' section in topology file of moleculetype 'molname'
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
    """Generate necessary strings for defining a BD3 type interaction in a topology

    :param int_dict: dictionary, 'BD3'  interaction dictionary for molname
    :param molname: molname
    :param bonded_tabpot: dictionary, bonded tabulated potentials.
    :param bonds: boolean, for QBB type, used to determine whether return 'bonds' string or 'angles' string
    :return: string, either bonds or angles
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
    """Generate 'dihedrals' section string for moleculetype 'molname'

    :param int_dict: dictionary, 'DIH' interaction dictionary
    :return: string for 'dihedrals' section in topology file of moleculetype 'molname'
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
    """Generate 'cmap' section string for moleculetype 'molname'

    :param int_dict: dictionary, 'CDI' interaction dictionary for molname
    :return: string for 'cmap' section in topology file of moleculetype 'molname'
    """
    cmap_string = ""
    for cmap_type, cmap_type_dict in int_dict.items():
        for atnum_list in cmap_type_dict.items():
            for at1, at2, at3, at4, at5 in atnum_list:
                cmap_string += f"{at1:8}{at2:8}{at3:8}{at4:8}{at5:8}\n"
    return cmap_string


def _gen_exclusions(int_dict):
    """Generate 'exclusions' section string for moleculetype 'molname'

    :param int_dict: dictionary, 'EXC' interaction dictionary for molname
    :return: string for 'exclusions' section in topology file of moleculetype 'molname'
    """
    exclusions_string = ""
    for atom_list in int_dict:
        for atom in atom_list:
            if int(atom) != 0:
                exclusions_string += f'{atom:5}'
        exclusions_string += '\n'
    return exclusions_string


def _find_moleculetypes(topology):
    """Finds all [ moleculetype ] sections, taking into account commented sections.

    :param topology: string, topology template file
    :return: list of [moleculetype] sections, in format [((begin_match, end_match), matching_string), ...]
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
    """Generates entire correctly formatted topology file

    :param topology: string, topology template file
    :param locations: list, locations of [moleculetype] interactions, i.e. output of _find_moleculetypes()
    :param topology_strings_dict: dictionary, contains {molname: {'.top_keyword' : section_string},...}, i.e. output of _gen_bonded_section_strings()
    :return: string, properly formatted topology file to be written out
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
    """Generates properly formatted [ moleculetype ] section for each moleculetype

    :param topology: string, template topology file
    :param bonded: dictionary, contains {'.top_keyword' : section_string},...}
    :param molname: molname
    :return: string, completed topology section for moleculetype 'molname'
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
    with open(f'{write_to}', 'w') as f:
        f.write(final_topology)
