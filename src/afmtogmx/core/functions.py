import numpy as np
from copy import deepcopy
import re

"""This module contains several valuable functions for the generation of input files for MD"""

total_bonded_added = 0  # counts total number of bonded interactions added


### Below begins the section defining pairwise interactions ###
def exp(param_list, r):
    """This function takes in a parameter list [P1 (A), P2 (alpha), 0, 0] and returns a potential and a force. \
     There is no modifying of parameters to change units, so it is up to the user to ensure proper units are being used\
      each time the funtion is called."""
    a, alpha = param_list[0], param_list[1]
    exponent = -alpha * r
    potential = a * np.exp(exponent)
    force = -potential * -alpha
    return potential, force


def srd(param_list, r):
    """This function takes in a parameter list [P1, P2, P3, 0] and returns a potential and a force. \
     There is no modifying of parameters, so it is up to the user to ensure proper units are being used \
     each time the funtion is called."""
    P1, P2, P3 = param_list[0], abs(param_list[1]), param_list[2]
    potential = P1 / (np.power(r, P2) + P3 ** P2)
    force_numerator = (P1 * P2 * np.power(r, P2))
    force_denom = r * (np.power((np.power(r, P2) + P3 ** P2), 2))
    force = np.divide(force_numerator, force_denom, out=np.zeros_like(force_denom), where=force_denom != 0)
    return potential, force


def shtr(param_list, r):
    """This function takes in a parameter list [P1, P2, P3, 0] and returns a potential and a force. \
     There is no modifying of parameters, so it is up to the user to ensure proper units are being used \
     each time the funtion is called."""
    P1, P2, P3 = param_list[0], param_list[1], param_list[2]
    temp_r = np.where(r <= P3, r, P3)
    potential = P1 * (np.power(temp_r, -P2) - 1 / (P3 ** P2) + P2 * (temp_r - P3) / (P3 ** (P2 + 1)))
    power = np.power(temp_r, (P2 + 1))
    force_scalar = P1 * P2
    force = -force_scalar * ((-1 / power) + 1 / P3 ** (P2 + 1))
    return potential, force


def powe(param_list, r):
    """This function takes in a parameter list [P1, P2, 0, 0] and returns a potential and a force. \
     There is no modifying of parameters, so it is up to the user to ensure proper units are being used \
     each time the function is called."""
    P1, P2, = param_list[0], param_list[1]
    potential = P1 * np.power(r, P2)
    force = -P1 * P2 * np.power(r, (P2 - 1))
    return potential, force

def quarbond(param_list, r):
    """This function takes in a parameter list [P1, P2, P3, P4] and returns a potential and a force. \
    There is no modifying of parameters, so it is up to the user to ensure proper units are being used \
    each time the function is called. The potential and force returned are numpy arrays."""
    P1, P2, P3, P4 = param_list[0], param_list[1], param_list[2], param_list[3]
    potential = (P2/2) * np.power(r - P1, 2) + (P3/3) * np.power(r - P1, 3) + (P4/4) * np.power(r - P1, 4)
    force = -(P2 * (r - P1) + P3 * np.power(r - P1, 2) + P4 * np.power(r - P1, 3))
    return potential, force



def gen_empty_bonded():
    """This function will generate an empty verion of the bonded dict necessary for the production of input files. This\
     empty bonded dict is populated later on in the workflow"""
    molname = dict()

    ### Start with bonded
    atoms = dict()
    atoms["All"] = dict()
    atoms["Virtual"] = dict()
    molname["ATO"] = atoms

    bonds = dict()
    bonds["HAR"] = dict()
    bonds["QUA"] = dict()
    molname["BON"] = bonds

    angles = dict()
    angles["HAR"] = dict()
    angles["QUA"] = dict()
    molname["ANG"] = angles

    bd3 = dict()
    bd3["QBB"] = dict()
    bd3["MUB"] = dict()
    molname["BD3"] = bd3

    dihedrals = dict()
    dihedrals["HAR"] = dict()
    dihedrals["NCO"] = dict()
    dihedrals["COS"] = dict()
    molname["DIH"] = dihedrals

    cdihedrals = dict()
    cdihedrals["CNCO"] = dict()
    cdihedrals["CCOS"] = dict()
    molname["CDI"] = cdihedrals

    exclusions = list()
    molname["EXC"] = exclusions

    return molname


def gen_empty_nonbonded():
    """This function will generate an empty verion of the nonbonded dict necessary for the production of input files.
    This empty nonbonded dict is populated later on in the workflow"""

    atom_pair = dict()

    nonbonded_int = ('COU', 'THC', 'GLJ', 'BUC', 'DBU', 'STR', 'EXP', 'POW', 'PEX', 'DPO', 'SRD')

    empty_list = []

    for item in nonbonded_int:
        atom_pair[f'{item}'] = empty_list

    return atom_pair


def _find_off_keywords(off_file_str=str):
    """Takes in an off file as a string and outputs a dictionary with each section.

    This script considers each off file to be made up of 5 sections: ff_input, intra_potential, inter_potential,
    molecular_definition, and table_potential. These 5 sections will be found and stored in a dictionary with the format
    {'section_name': 'section'}. This sections dictionary will be returned to the user after successful execution of this
    script.

    :param off_file_str: str
    :return: dictionary containing section names, and sections for each of the 5 sections in an .off file
    """
    ff_input = off_file_str.find("Atom Types:")
    intra_potential = off_file_str.find("Intra-Potential:\n")
    inter_potential = off_file_str.find("Inter-Potential:\n")
    molecular_definition = off_file_str.find("Molecular-Definition:\n")
    table_potential = off_file_str.find("Table-Potential:\n")
    file_end = len(off_file_str) - 1
    sections = {'ff_input': (0, ff_input),
                'intra_potential': (intra_potential, inter_potential),
                'inter_potential': (inter_potential, molecular_definition),
                'molecular_definition': (molecular_definition, table_potential),
                'table_potential': (table_potential, file_end)
                }
    for key, value in sections.items():
        sections[key] = off_file_str[value[0]:value[1]]

    return sections


def _recognize_keywords(section=str):
    """Takes in a section str and outputs a list containing [['key', beginning, ending]...]

    :param section: string of off file section
    :return: list of lists
    """
    import re
    total_keywords = ('OPT', 'KEY', 'MOL', 'ATO', 'BON', 'ANG', 'BD3', 'DIH', 'CDIH', 'EXC', 'FUD', 'COU', 'THC', 'GLJ',
                      'BUC', 'DBU', 'STR', 'EXP', 'POW', 'PEX', 'DPO', 'SRD', 'EQV', 'CHA')

    test = re.findall('\[.*]', section)
    itert = re.finditer('\[.*]', section)
    keywords_and_locations = []
    for item in itert:
        letters = re.findall(r'[a-zA-Z0-9]', item[0])
        three_str = ''.join(letters[:3])
        four_str = ''.join(letters[:4])
        if three_str in total_keywords:
            keywords_and_locations.append([three_str, item.start(), item.end()])
            three_str = None
        elif four_str in total_keywords:
            keywords_and_locations.append([four_str, item.start(), item.end()])
            four_str = None

    return keywords_and_locations


def _filter_interactions(interactions=list):
    """Filters list containing keywords and locations in sections[ff_input]. Ignores all other keywords

    :param interactions: list of keywords and locations
    :return: bonded_list, nonbonded_list
    """
    bonded_keywords = ('MOL', 'ATO', 'BON', 'ANG', 'BD3', 'DIH', 'CDIH', 'EXC')
    nonbonded_keywords = ('FUD', 'COU', 'THC', 'GLJ', 'BUC', 'DBU', 'STR', 'EXP', 'POW', 'PEX', 'DPO', 'SRD', 'EQV',
                          'CHA')
    bonded_list = []
    nonbonded_list = []
    for item in interactions:
        if item[0] in bonded_keywords:
            bonded_list.append(item)
        elif item[0] in nonbonded_keywords:
            nonbonded_list.append(item)

    return bonded_list, nonbonded_list


def _find_end_bonded(final_bonded=list, ff_input=str):
    """Finds appropriate location to end reading the bonded section of ff_input

    :param final_bonded: last bonded term in list containing keywords and locations of bonded interactions
    :param ff_input: top of off file, original input file to cryoff
    :return: list containing ['END', end_location, end_location]
    """
    temp_str = ff_input[final_bonded[-1]:]
    end_location = re.search('\[.*]', temp_str)
    end_location = end_location.start() + final_bonded[-1]
    return ['END', end_location, end_location]


def _find_molnames(bonded=list, ff_input=str):
    """Finds molnames based on keywords and locations of bonded interactions, and ff_input

    :param bonded: list of keywords and locations of all bonded interactions
    :param ff_input: top of off file
    :return: list of molnames
    """
    molnames = []
    for item in bonded:
        if item[0] != 'MOL':
            pass
        else:
            begin = item[2]
            end = re.search('\n', ff_input[item[2]:]).start() + begin
            molname = ff_input[begin:end].strip()
            molname = molname.split()[0]
            molnames.append(molname)

    return molnames


def _split_into_molecules(bonded=list):
    """Takes in a list of all bonded interactions and splits into list of lists where each individual sublist has 1 MOL.

    :param bonded: list containing all the bonded terms in sections['ff_input']
    :return: list of lists; each sublist corresponds to 1 molecule
    """
    result = []
    sublist = []
    for item in bonded:
        if item[0] == 'MOL':
            if sublist:  # Append the sublist if it's not empty
                sublist.append(['END', item[1], item[1]])
                result.append(sublist)
                sublist = [item]  # Start a new sublist with the current 'MOL'
            else:
                sublist.append(item)  # Include 'MOL' in an empty sublist
        else:
            sublist.append(item)
    if sublist:  # Append the last sublist if it's not empty
        result.append(sublist)

    return result


def _gather_fitted_bonded(ff_input=str):  # splits sections['intra_potential'] into list of strings based on each
    # different [ TERM ] section
    return re.split(r'\n[ ]+\[.+]', ff_input)[1:]


def _parse_bonded(unsorted_bonded, molecule=list, ff_input=str):
    """Begins work of actually splitting the molecules and sorting through which parameters go where

    :param unsorted_bonded: list of sections['intra_potential'] split up by each [ INTERACTION ]
    :param molecule: list containing all sections related to 1 MOL, with [['INTERACTION', beginning, ending]...,['END'..
    :param ff_input: sections['ff_input'], a.k.a. top of .off file corresponding to .ff file
    :return: dictionary with all elements related to bonded portion of each specified molecule
    """
    populated_molecule = gen_empty_bonded()  # initialize empty dictionary with proper structure
    for num in range(0, len(molecule)):  # iterate over molecule
        key = molecule[num][0]
        if key == 'END':  # If the key is 'END', stop iterating
            break
        beg = molecule[num][2]  # ending of current section to
        end = molecule[num + 1][1]  # beginning of next section, gives the true section information as a string

        if key == 'MOL':  # Do nothing if 'MOL', already extracted information from here
            continue

        # call function to do heavy lifting of breaking each key into proper sections and saving to dictionary
        populated_molecule[f'{key}'] = _parse_bonded_section(unsorted_bonded, section=key, beginning=beg, ending=end,
                                                             ff_input=ff_input)
    return populated_molecule


def _parse_bonded_section(uns_bonded, section=str, beginning=int, ending=int, ff_input=str):
    """Does the heavy lifting for sorting through and properly formatting dictionaries for bonded terms.

    :param uns_bonded: output of _gather_fitted_bonded
    :param section: which term to parse
    :param beginning: beginning of section
    :param ending: end of section
    :param ff_input: sections['ff_input'], a.k.a .ff portion at top of .off file
    :return: portion of dictionary for each 'key', properly formatted with fitted parameters
    """
    section_str = ff_input[beginning:ending]
    global total_bonded_added

    if section == "ATO":
        All = dict()  # Holds every atom
        Virtual = dict()  # Holds only virtual atoms
        lines = section_str.split("\n")[1:-1]  # split lines up, keeping only non-blank lines
        for line in lines:
            if "*" in line:  # if virtual atom
                nostar = re.sub("\*", "", line)
                atnum, at1, at2 = nostar.split()[:3]
                atnum = int(atnum)
                definition = tuple(nostar.split()[3:])
                Virtual[(atnum, at1, at2)] = definition  # save definition of virtual atom location for later just
                # in case it will be needed
                All[atnum] = (at1, at2)  # add atnum, at1, at2 to all
            else:
                atnum, at1, at2 = line.split()[:3]
                All[atnum] = (at1, at2)  # add atnum, at1, at2 to all
        ATO_dict = {'All': All, 'Virtual': Virtual}
        return ATO_dict  # properly formatted 'ATO' sub-dictionary

    if section == "BON":
        HAR = dict()
        QUA = dict()

        lines = section_str.split("\n")[1:-1]
        is_qua = False
        is_har = False
        for line in lines:
            if "QUA" in line:
                total_bonded_added += 1  # iterate counter to keep track of which parameters to use from uns_bonded
                is_qua, is_har = True, False
                P1, P2, P3, P4 = [float(i) for i in uns_bonded[total_bonded_added - 1].split("\n")[0].split()[-4:]]
                QUA[(P1, P2, P3, P4)] = []
                continue
            if "HAR" in line:
                total_bonded_added += 1
                is_har, is_qua = True, False
                P1, P2, = [float(i) for i in uns_bonded[total_bonded_added - 1].split("\n")[0].split()[-2:]]
                HAR[(P1, P2)] = []
                continue
            at1, at2 = [int(i) for i in line.split()]  # split up atoms
            if is_qua:
                QUA[(P1, P2, P3, P4)].append([at1, at2])  # add to QUA
            if is_har:
                HAR[(P1, P2)].append([at1, at2])  # add to HAR
        BON_dict = {'HAR': HAR, 'QUA': QUA}  # final, properly formatted 'BON' sub-dictionary

        return BON_dict

    if section == "ANG":  # Essentially identical to BON, see that section
        HAR = dict()
        QUA = dict()
        lines = section_str.split("\n")[1:-1]
        is_qua = False
        is_har = False
        for line in lines:
            if "QUA" in line:
                total_bonded_added += 1
                is_qua, is_har = True, False
                P1, P2, P3, P4 = [float(i) for i in uns_bonded[total_bonded_added - 1].split("\n")[0].split()[-4:]]
                QUA[(P1, P2, P3, P4)] = []
                continue
            if "HAR" in line:
                total_bonded_added += 1
                is_har, is_qua = True, False
                P1, P2 = [float(i) for i in uns_bonded[total_bonded_added - 1].split("\n")[0].split()[-2:]]
                HAR[(P1, P2)] = []
                continue
            at1, at2, at3 = [int(i) for i in line.split()]
            if is_qua:
                QUA[(P1, P2, P3, P4)].append([at1, at2, at3])
            if is_har:
                HAR[(P1, P2)].append([at1, at2, at3])

        ANG_dict = {'HAR': HAR, 'QUA': QUA}

        return ANG_dict

    if section == "BD3":  # NEEDS TESTING, JUST GUESSING AT THE MOMENT
        QBB = dict()
        MUB = dict()
        lines = section_str.split("\n")[1:-1]
        is_qbb = False
        is_mub = False
        for line in lines:
            if "QBB" in line:
                total_bonded_added += 1
                is_qbb, is_mub = True, False
                P1, P2, P3, P4, P5 = [float(i) for i in uns_bonded[total_bonded_added - 1].split("\n")[0].split()[-5:]]
                QBB[(P1, P2, P3, P4, P5)] = []
                continue
            if "MUB" in line:
                total_bonded_added += 1
                is_mub, is_qbb = True, False
                P1, P2, P3, P4 = [float(i) for i in uns_bonded[total_bonded_added - 1].split("\n")[0].split()[-4:]]
                continue
            at1, at2, at3 = [int(i) for i in line.split()]
            if is_qbb:
                QBB[(P1, P2, P3, P4, P5)].append([at1, at2, at3])
            if is_mub:
                MUB[(P1, P2, P3, P4)].append([at1, at2, at3])
        BD3_dict = {'QBB': QBB, 'MUB': MUB}

        return BD3_dict

    if section == "DIH":
        HAR = dict()
        NCO = dict()
        COS = dict()
        lines = section_str.split("\n")[1:-1]
        is_har = False
        is_nco = False
        is_cos = False
        for line in lines:
            if "HAR" in line:
                total_bonded_added += 1
                is_har, is_nco, is_cos = True, False, False
                P1, P2 = [float(i) for i in uns_bonded[total_bonded_added - 1].split("\n")[0].split()[-2:]]
                HAR[(P1, P2)] = []
                continue

            if "NCO" in line:
                total_bonded_added += 1
                is_nco, is_har, is_cos = True, False, False
                P1, P2, P3 = [float(i) for i in uns_bonded[total_bonded_added - 1].split("\n")[0].split()[-3:]]
                NCO[(P1, P2, P3)] = []
                continue

            if "COS" in line:
                total_bonded_added += 1
                is_cos, is_har, is_nco = True, False, False
                P1, P2, P3 = [float(i) for i in uns_bonded[total_bonded_added - 1].split("\n")[0].split()[-3:]]
                COS[(P1, P2, P3)] = []
                continue

            at1, at2, at3, at4 = [int(i) for i in line.split()]
            if is_har:
                HAR[(P1, P2)].append([at1, at2, at3, at4])
            if is_nco:
                NCO[(P1, P2, P3)].append([at1, at2, at3, at4])
            if is_cos:
                COS[(P1, P2, P3)].append([at1, at2, at3, at4])

        DIH_dict = {'HAR': HAR, 'NCO': NCO, 'COS': COS}

        return DIH_dict

    if section == "CDIH":  # NEEDS TESTING
        CNCO = dict()
        CCOS = dict()
        lines = section_str.split("\n")[1:-1]
        is_cnco = False
        is_ccos = False
        for line in lines:
            if "CNCO" in line:
                total_bonded_added += 1
                is_cnco, is_ccos = True, False
                P1, P2, P3, P4 = [float(i) for i in uns_bonded[total_bonded_added - 1].split("\n")[0].split()[-4:]]
                CNCO[(P1, P2, P3, P4)] = []
                continue
            if "CCOS" in line:
                total_bonded_added += 1
                is_ccos, is_cnco = True, False
                P1, P2, P3, P4 = [float(i) for i in uns_bonded[total_bonded_added - 1].split("\n")[0].split()[-4:]]
                CCOS[(P1, P2, P3, P4)] = []
                continue

            at1, at2, at3, at4 = [int(i) for i in line.split()]
            if is_cnco:
                CNCO[(P1, P2, P3, P4)].append([at1, at2, at3, at4])
            if is_ccos:
                CCOS[(P1, P2, P3, P4)].append([at1, at2, at3, at4])

        CDIH_dict = {"CNCO": CNCO, "CCOS": CCOS}

        return CDIH_dict

    if section == "EXC":
        EXC = list()
        lines = section_str.split("\n")[1:-1]
        for line in lines:
            EXC.append([int(i) for i in line.split()])
        return EXC


def _clean_inter_potential(inter_potential=str):
    """Splits up inter potential into list of lists with format: [[('At1', 'At2'), 'Interaction', P1, P2,...]...]

    :param inter_potential:
    :return: list of lists containing pair, interaction, and parameters for each item in inter_potential section
    """
    inter_potential = re.sub('\:', '', inter_potential)  # remove colon from inter_potential
    inter_potential = inter_potential.split("\n")[1:-2]  # split inter potential by lines
    inter_potential = [i.split()[:-4] for i in inter_potential]  # split inter potential,
    # keeping only At1~At1, INT, Parameters
    # Organize according to predetermined rule
    for item in inter_potential:
        item[0] = tuple(sorted(item[0].split("~")))
        item[1] = item[1]
#        item[1] = item[1][:3]

        for num in range(2, len(item)):
            item[num] = float(item[num])
    return inter_potential

def _gen_empty_charge_dict(bonded):
    """Generates empty charge dictionary with format {"Mol1" : {'At1' : 0.0, 'At2' : 0.0...}, 'Mol2':...}

    :param bonded: self.bonded dictionary
    :return: dictionary with molnames atoms and charges, with all charges equal to 0.0
    """
    charges = {}
    for mol in bonded:
        all_atoms = {pair[1]: 0.0 for key, pair in bonded[mol]['ATO']['All'].items() if
                     pair[1] != 'NETF' and pair[1] != 'TORQ'}  # setting each atom in each mol charge == 0.0,placeholder
        charges[mol] = all_atoms  # set charges[mol] equal to what we just defined and return
    return charges


def _gen_charges_from_known(charges, nonbonded, known_atom, known_atom_charge):
    """Generates charges for each atom based on a known atom charge, often times

    :param charges: dictionary containing molname, atoms, and 0.0 as charges
    :param nonbonded: self.nonbonded dictionary
    :param known_atom: known_atom, str
    :param known_atom_charge: charge, float
    :return: populated charges dictionary
    """
    temp_charges = {}
    temp_charges[known_atom] = known_atom_charge  # set known_atom_charge

    for pair, value in nonbonded.items():
        # iterate over pairs, and charges, iff known_atom in pair and not (known, known)
        if 'COU' in value and pair != (known_atom, known_atom) and known_atom in pair:
            list_pair = list(pair)
            list_pair.remove(known_atom)  # remove known atom from pair, calculate charge
            charge = value['COU'][0][0] / known_atom_charge
            temp_charges[list_pair[0]] = charge  # set temp_charges[Atom] to charge

    for molname, atoms in charges.items():
        for atom in atoms:
            if atom in temp_charges:
                charges[molname][atom] = temp_charges[atom]  # assign charges to each atom in charges, per molname
    return charges


def _get_known_atom_charge(known_atom, nonbonded, sign):
    """If no atom charge is known, calculate charges based on single ('known_atom', 'known_atom') pair, with some
    nonzero, positive or negative charge, via (+/-)sqrt(('known_atom', 'known_atom'))

    :param known_atom: atom with a nonzero coulombic interaction with its self
    :param nonbonded: self.nonbonded dictionary
    :param sign: sign of final individual atom, either '+' or '-'
    :return: derived charge of single selected known_atom
    """
    try:  # try to get pair charge
        pair_charge = nonbonded[(known_atom, known_atom)]['COU'][0][0]
    except KeyError as err:  # unless there are no (known_atom, known_atom) with COUlomb
        print("Must use a known atom that has an Atom ~ Atom COUlomb interaction in the .off file")
        exit(1)
    if pair_charge == 0.0:  # If pair charge = 0, can't use
        print("Must use a known atom that has a nonzero charge in the Atom ~ Atom COUlomb interaction in the .off file")
        exit(1)
    if sign == '+':  # set single_charge based on sign
        print(pair_charge)
        single_charge = pair_charge ** (1 / 2)
    elif sign == '-':
        single_charge = -pair_charge ** (1 / 2)
    else:
        print("Must use '+' or '-' when stating sign value for known_atom")
        exit(1)

    return single_charge


def _normalize_charges(normalization, charges, bonded, tolerance):
    """Normalizes charges if it is found that leftover charge is greater than some tolerance value

    :param normalization: normalization method. Current options are only 'M-POPULOUS'
    :param charges: charges dictionary with molname, atoms, charges per atom
    :param bonded: self.bonded
    :param tolerance: default of 1E-5
    :return: charges dictionary with normalized, neutral molecules
    """
    updated_charges = {}
    for mol in bonded:
        total_charge = 0
        all_atoms = [pair[1] for key, pair in bonded[mol]['ATO']['All'].items() if
                     pair[1] != 'NETF' and pair[1] != 'TORQ']  # holds every atom per molname
        for atom in all_atoms:
            total_charge += charges[mol][atom]  # calculates total leftover charge
        if abs(total_charge) > tolerance:
            print(f"Modifying charges for molecule {mol}")
            print(f"Initial Total Charge: {total_charge}")
            # pass to function _normalization_process to perform heavy lifting
            updated_charges[mol] = _normalization_process(normalization, all_atoms, charges[mol], total_charge)
        else:
            updated_charges[mol] = charges[mol]  # if already within tolerance, keep the same charges

    return updated_charges


def _normalization_process(normalization, all_atoms, charges, total_charge):
    """Performs heavy lifting of normalizing/neutralizing molecules if above the tolerance level

    :param normalization: Normalization method; currently only support M-POPULOUS. M-POPULOUS divides the leftover
    charge by the number of most populous atoms, and then subtracts leftover_charge/num(most populous) from the most
    populous atom charge, to give effectively a zero charge for the molecule.
    :param all_atoms: list of all atoms, with repeats, from each molname
    :param charges: charges dictionary
    :param total_charge: total charge before modifying charges per molecule
    :return: charges dictionary per molname with neutralized molecule
    """

    if normalization == "M-POPULOUS":
        from collections import defaultdict
        atom_instances = defaultdict(int)  # count number of atoms per molecule
        abs_total_charge = 0
        for atom in all_atoms:
            atom_instances[atom] += 1
            abs_total_charge += abs(charges[atom])

        atom_instances = dict(sorted(atom_instances.items(), key=lambda item: item[1], reverse=True))  # sort number of
        # atoms per molecule as from highest number to lowest number of atoms

        for atom in atom_instances:
            if charges[atom] != 0:  # when nonzero charge is found
                num_atoms = atom_instances[atom]
                subtract_charge = total_charge / num_atoms
                charges[atom] -= subtract_charge  # subtract leftover charge/num atoms from original atom charge
                print(f"Modifying atom {atom} by {-subtract_charge}")
                break
        new_charge = 0
        for atom in all_atoms:
            new_charge += charges[atom]  # calculate the new leftover charge, will be effectively zero
        print("New Charge: ", new_charge)
        print("")
        afternorm_charge = 0
        for atom in all_atoms:
            afternorm_charge += round(charges[atom], 5)
        print(afternorm_charge)
        return charges  # return neutralized dictionary

def _filter_bonded(bonded):
    to_filter = ('BON', 'ANG', 'BD3', 'DIH', 'CDI', 'EXC')
    filtered = dict()

    for interaction, int_dict in bonded.items():
        if interaction == 'EXC':  # skip EXC
            filtered[interaction] = deepcopy(int_dict)
        else:  # Below is section of code to quickly remove empty portions of bonded. Also removes 'Virtual' if no virtual atoms
            for bonded_term in int_dict:
                if int_dict[bonded_term]:
                    if interaction in filtered:
                        filtered[interaction][bonded_term] = deepcopy(int_dict[bonded_term])
                    else:
                        filtered[interaction] = dict()
                        filtered[interaction][bonded_term] = deepcopy(int_dict[bonded_term])
    return filtered


def _remove_empty_and_cou_interactions_nonbonded(nonbonded):
    """Takes in self.nonbonded and returns an identical dictionary with COU and empty interactions removed

    :param nonbonded: self.nonbonded
    :return: self.nonbonded with empty interactions remove
    """
    cleaned_dict = {}
    for pair, int_dict in nonbonded.items():
        pair_dict = {k: v for k, v, in int_dict.items() if v}
        if 'COU' in pair_dict:
            del pair_dict['COU']
        if pair_dict:
            cleaned_dict[pair] = pair_dict
    return cleaned_dict