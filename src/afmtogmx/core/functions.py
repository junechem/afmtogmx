import numpy as np
from copy import deepcopy
import re
from collections import defaultdict

"""This module contains several valuable functions for the generation of input files for MD"""

total_bonded_added = 0  # counts total number of bonded interactions added


### Below begins the section defining pairwise interactions ###
def exp(param_list, r):
    """Calculate the potential and force of an exponential interaction.

    Parameters
    ----------
    param_list : list or tuple
        A list of parameters for the exponential function, expected to be
        in the format `[A, alpha, 0, 0]`.
    r : float or numpy.ndarray
        The distance(s) at which to calculate the potential and force.

    Returns
    -------
    tuple of (float or numpy.ndarray, float or numpy.ndarray)
        A tuple containing the calculated potential and force.

    Notes
    -----
    - The function calculates the potential as `U(r) = A * exp(-alpha * r)`.
    - The force is calculated as the negative derivative of the potential.
    - No unit conversions are performed within this function. The user is
      responsible for ensuring that the input parameters and distances
      have consistent units.
    """
    a, alpha = param_list[0], param_list[1]
    exponent = -alpha * r
    potential = a * np.exp(exponent)
    force = -potential * -alpha
    return potential, force


def srd(param_list, r):
    """Calculate the potential and force of a Short-Range Damped (SRD) interaction.

    The SRD potential is defined as :math:`U(r) = P_1 / (r^{P_2} + P_3^{P_2})`.

    Parameters
    ----------
    param_list : list or tuple
        A list of parameters for the SRD function, expected to be
        in the format `[P1, P2, P3, 0]`. Note that P2 is taken as its absolute value.
    r : float or numpy.ndarray
        The distance(s) at which to calculate the potential and force.

    Returns
    -------
    tuple of (float or numpy.ndarray, float or numpy.ndarray)
        A tuple containing the calculated potential and force.

    Notes
    -----
    - The force is calculated as the negative derivative of the potential.
    - No unit conversions are performed within this function. The user is
      responsible for ensuring that the input parameters and distances
      have consistent units.
    - Handles potential division by zero for force calculation by
      returning zero where `force_denom` is zero.
    """
    P1, P2, P3 = param_list[0], abs(param_list[1]), param_list[2]
    potential = P1 / (np.power(r, P2) + P3 ** P2)
    force_numerator = (P1 * P2 * np.power(r, P2))
    force_denom = r * (np.power((np.power(r, P2) + P3 ** P2), 2))
    force = np.divide(force_numerator, force_denom, out=np.zeros_like(force_denom), where=force_denom != 0)
    return potential, force


def shtr(param_list, r):
    """Calculate the potential and force of a Shifted Truncated (SHTR) interaction.

    Parameters
    ----------
    param_list : list or tuple
        A list of parameters for the SHTR function, expected to be
        in the format `[P1, P2, P3, 0]`.
    r : float or numpy.ndarray
        The distance(s) at which to calculate the potential and force.

    Returns
    -------
    tuple of (float or numpy.ndarray, float or numpy.ndarray)
        A tuple containing the calculated potential and force.

    Notes
    -----
    - The potential is defined for `r <= P3` and truncated/shifted such
      that `U(P3) = 0` and `U'(P3) = 0`.
    - The calculation uses `temp_r = np.where(r <= P3, r, P3)`.
    - No unit conversions are performed within this function. The user is
      responsible for ensuring that the input parameters and distances
      have consistent units.
    """
    P1, P2, P3 = param_list[0], param_list[1], param_list[2]
    temp_r = np.where(r <= P3, r, P3)
    potential = P1 * (np.power(temp_r, -P2) - 1 / (P3 ** P2) + P2 * (temp_r - P3) / (P3 ** (P2 + 1)))
    power = np.power(temp_r, (P2 + 1))
    force_scalar = P1 * P2
    force = -force_scalar * ((-1 / power) + 1 / P3 ** (P2 + 1))
    return potential, force


def powe(param_list, r):
    """Calculate the potential and force of a power-law interaction.

    The power-law potential is defined as :math:`U(r) = P_1 * r^{P_2}`.

    Parameters
    ----------
    param_list : list or tuple
        A list of parameters for the power-law function, expected to be
        in the format `[P1, P2, 0, 0]`.
    r : float or numpy.ndarray
        The distance(s) at which to calculate the potential and force.

    Returns
    -------
    tuple of (float or numpy.ndarray, float or numpy.ndarray)
        A tuple containing the calculated potential and force.

    Notes
    -----
    - The force is calculated as the negative derivative of the potential.
    - No unit conversions are performed within this function. The user is
      responsible for ensuring that the input parameters and distances
      have consistent units.
    """
    P1, P2, = param_list[0], param_list[1]
    potential = P1 * np.power(r, P2)
    force = -P1 * P2 * np.power(r, (P2 - 1))
    return potential, force

def quarbond(param_list, r):
    """Calculate the potential and force for a quartic bond interaction.

    The quartic bond potential is defined as:
    :math:`U(r) = (P_2/2)(r - P_1)^2 + (P_3/3)(r - P_1)^3 + (P_4/4)(r - P_1)^4`

    Parameters
    ----------
    param_list : list or tuple
        A list of parameters for the quartic bond function, expected to be
        in the format `[P1, P2, P3, P4]`.
    r : float or numpy.ndarray
        The distance(s) at which to calculate the potential and force.

    Returns
    -------
    tuple of (numpy.ndarray, numpy.ndarray)
        A tuple containing the calculated potential and force, both as NumPy arrays.

    Notes
    -----
    - The force is calculated as the negative derivative of the potential.
    - No unit conversions are performed within this function. The user is
      responsible for ensuring that the input parameters and distances
      have consistent units.
    """
    P1, P2, P3, P4 = param_list[0], param_list[1], param_list[2], param_list[3]
    potential = (P2/2) * np.power(r - P1, 2) + (P3/3) * np.power(r - P1, 3) + (P4/4) * np.power(r - P1, 4)
    force = -(P2 * (r - P1) + P3 * np.power(r - P1, 2) + P4 * np.power(r - P1, 3))
    return potential, force



def gen_empty_bonded():
    """Generate an empty dictionary structure for bonded interactions.

    This function creates a nested dictionary that serves as a template
    for storing all bonded interaction parameters parsed from an .off file.
    This empty structure is designed to be populated later in the workflow.

    Returns
    -------
    dict
        An empty dictionary with the following structure, representing a
        single molecule's bonded interactions:
        {
            "ATO": {"All": {}, "Virtual": {}},
            "BON": {"HAR": {}, "QUA": {}},
            "ANG": {"HAR": {}, "QUA": {}},
            "BD3": {"QBB": {}, "MUB": {}},
            "DIH": {"HAR": {}, "NCO": {}, "COS": {}},
            "CDI": {"CNCO": {}, "CCOS": {}},
            "EXC": []
        }
    """
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
    """Generate an empty dictionary structure for nonbonded interactions.

    This function creates a dictionary that serves as a template for
    storing nonbonded interaction parameters for a given atom pair.
    The empty structure is designed to be populated later in the workflow.

    Returns
    -------
    dict
        An empty dictionary with keys representing various nonbonded
        interaction types (e.g., 'COU', 'THC', 'GLJ', 'BUC', etc.),
        each mapped to an empty list.
        Example: `{'COU': [], 'THC': [], 'GLJ': [], ...}`
    """

    atom_pair = dict()

    nonbonded_int = ('COU', 'THC', 'GLJ', 'BUC', 'DBU', 'STR', 'EXP', 'POW', 'PEX', 'DPO', 'SRD')

    empty_list = []

    for item in nonbonded_int:
        atom_pair[f'{item}'] = empty_list

    return atom_pair


def _find_off_keywords(off_file_str=str):
    """Parse an .off file string into its constituent sections.

    This function takes the content of an .off file as a string and
    identifies five predefined sections: `ff_input`, `intra_potential`,
    `inter_potential`, `molecular_definition`, and `table_potential`.
    Each section's content is extracted and stored in a dictionary.

    Parameters
    ----------
    off_file_str : str
        The entire content of an .off file as a single string.

    Returns
    -------
    dict
        A dictionary where keys are the section names (e.g., 'ff_input')
        and values are the string content of each corresponding section.
        Example: `{'ff_input': 'Atom Types:\\n...', 'intra_potential': 'Intra-Potential:\\n...', ...}`
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
    """Identify and locate keywords within a given section string.

    This function scans a section of an .off file (provided as a string)
    to find predefined keywords and their starting and ending character
    locations.

    Parameters
    ----------
    section : str
        A string representing a section of an .off file.

    Returns
    -------
    list of list
        A list where each inner list contains `[keyword, start_index, end_index]`.
        `keyword` is a string matching one of the predefined keywords.
        `start_index` and `end_index` are integers indicating the keyword's
        span within the `section` string.

    Notes
    -----
    The function searches for keywords including:
    ('OPT', 'KEY', 'MOL', 'ATO', 'BON', 'ANG', 'BD3', 'DIH', 'CDIH', 'EXC',
    'FUD', 'COU', 'THC', 'GLJ', 'BUC', 'DBU', 'STR', 'EXP', 'POW', 'PEX',
    'DPO', 'SRD', 'EQV', 'CHA')
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
    """Filter a list of interactions into bonded and nonbonded categories.

    This function takes a list of identified keywords and their locations
    from the `ff_input` section of an .off file and separates them into
    two distinct lists: one for bonded interactions and one for nonbonded
    interactions, based on predefined keyword sets.

    Parameters
    ----------
    interactions : list of list
        A list where each inner list contains `[keyword, start_index, end_index]`
        as identified by `_recognize_keywords`.

    Returns
    -------
    tuple of (list, list)
        A tuple containing two lists:
        - `bonded_list`: A list of interactions identified as bonded.
        - `nonbonded_list`: A list of interactions identified as nonbonded.

    Notes
    -----
    - **Bonded Keywords**: ('MOL', 'ATO', 'BON', 'ANG', 'BD3', 'DIH', 'CDIH', 'EXC')
    - **Nonbonded Keywords**: ('FUD', 'COU', 'THC', 'GLJ', 'BUC', 'DBU', 'STR',
                              'EXP', 'POW', 'PEX', 'DPO', 'SRD', 'EQV', 'CHA')
    Any other keywords found in the input `interactions` list will be ignored.
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
    """Determine the end location for reading the bonded section of the force field input.

    This function helps to correctly delimit the bonded section within the
    force field input string by finding the position of the next keyword
    after the last identified bonded interaction.

    Parameters
    ----------
    final_bonded : list
        A list containing the last bonded interaction term, typically
        in the format `[keyword, start_index, end_index]`.
    ff_input : str
        The string content of the `ff_input` section of the .off file.

    Returns
    -------
    list
        A list in the format `['END', end_location, end_location]`, where
        `end_location` is the character index marking the end of the
        bonded section.
    """
    temp_str = ff_input[final_bonded[-1]:]
    end_location = re.search('\[.*]', temp_str)
    end_location = end_location.start() + final_bonded[-1]
    return ['END', end_location, end_location]


def _find_molnames(bonded=list, ff_input=str):
    """Extract molecule names from the `ff_input` section of an .off file.

    This function identifies molecule definitions (`MOL` keywords) within
    the bonded interactions list and extracts the corresponding molecule
    names from the `ff_input` string.

    Parameters
    ----------
    bonded : list of list
        A list containing keyword-location information for all bonded
        interactions, as generated by `_filter_interactions`.
    ff_input : str
        The string content of the `ff_input` section of the .off file.

    Returns
    -------
    list of str
        A list of extracted molecule names.
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
    """Split a list of all bonded interactions into molecule-specific sublists.

    This function processes a comprehensive list of all bonded interactions,
    identified by `MOL` keywords, and organizes them into a list of lists.
    Each inner list represents the bonded interactions pertinent to a
    single molecule. An 'END' marker is inserted to delimit each molecule's
    interactions.

    Parameters
    ----------
    bonded : list of list
        A list containing all bonded terms from the `ff_input` section,
        including `MOL` keywords.

    Returns
    -------
    list of list
        A list where each sublist corresponds to a single molecule and
        contains its specific bonded interaction terms.
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


def _gather_fitted_bonded(ff_input=str):
    """Split the `intra_potential` section into individual bonded interaction strings.

    This function takes the string content of the `intra_potential` section
    from an .off file and divides it into a list of strings, where each
    string corresponds to a different bonded interaction defined within that section.

    Parameters
    ----------
    ff_input : str
        The string content of the `intra_potential` section of the .off file.

    Returns
    -------
    list of str
        A list of strings, where each string represents a fitted bonded
        interaction block (e.g., a `[ TERM ]` section) from the
        `intra_potential`.
    """    return re.split(r'\n[ ]+\[.+]', ff_input)[1:]


def _parse_bonded(unsorted_bonded, molecule=list, ff_input=str):
    """Parse bonded interactions for a single molecule.

    This function takes a list of a molecule's bonded interactions and
    the overall `intra_potential` section to extract and organize the
    fitted parameters for that specific molecule into a structured dictionary.
    It uses `gen_empty_bonded` to initialize the dictionary and
    `_parse_bonded_section` for detailed parsing.

    Parameters
    ----------
    unsorted_bonded : list of str
        A list of strings, where each string represents a fitted bonded
        interaction block from the `intra_potential` section, as output
        by `_gather_fitted_bonded`.
    molecule : list of list
        A list containing keyword-location information for all sections
        related to a single molecule, including an 'END' marker.
    ff_input : str
        The string content of the `ff_input` section of the .off file.

    Returns
    -------
    dict
        A dictionary containing all elements related to the bonded portion
        of the specified molecule, structured according to `gen_empty_bonded`.
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
    """Parse a specific bonded interaction section and format its dictionary.

    This function performs the detailed parsing of a single bonded interaction
    section (e.g., 'ATO', 'BON', 'ANG') from the `ff_input` based on its
    start and end locations. It extracts parameters and atom definitions,
    then populates and returns a structured dictionary for that specific
    interaction type.

    Parameters
    ----------
    uns_bonded : list of str
        A list of strings, where each string represents a fitted bonded
        interaction block from the `intra_potential` section, as output
        by `_gather_fitted_bonded`. This is used to retrieve parameter values.
    section : str
        The keyword string identifying the specific bonded interaction
        section to parse (e.g., 'ATO', 'BON', 'ANG').
    beginning : int
        The starting character index of the section within `ff_input`.
    ending : int
        The ending character index of the section within `ff_input`.
    ff_input : str
        The string content of the `ff_input` section of the .off file.

    Returns
    -------
    dict
        A dictionary containing the parsed and formatted parameters
        and atom lists for the specified `section`. The structure varies
        by section type (e.g., 'ATO' has 'All' and 'Virtual' keys,
        'BON' has 'HAR' and 'QUA' keys).
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
    """Parse and clean the `inter_potential` section into a structured list.

    This function takes the raw string content of the `inter_potential`
    section from an .off file, performs several cleaning and parsing steps,
    and returns a list of lists. Each inner list represents a nonbonded
    interaction with its atom pair, interaction type, and parameters.

    Parameters
    ----------
    inter_potential : str
        The string content of the `inter_potential` section of the .off file.

    Returns
    -------
    list of list
        A list where each inner list has the format:
        `[(At1, At2), 'InteractionType', P1, P2, ...]`.
        - `(At1, At2)`: A tuple of sorted atom names (strings) for the pair.
        - `InteractionType`: A string identifying the interaction.
        - `P1, P2, ...`: Floats representing the parameters for the interaction.

    Notes
    -----
    The processing steps include:
    - Removing colons (`:`) from the string.
    - Splitting the string into lines and then into components.
    - Sorting atom pairs alphabetically (e.g., `C~H` becomes `('C', 'H')`).
    - Converting interaction parameters to float types.
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


def _filter_bonded(bonded):
    """Remove empty interaction definitions from a bonded interactions dictionary.

    This function processes a dictionary of bonded interactions for a molecule,
    removing any interaction types or terms that are empty. It also handles
    the special case of 'EXC' (exclusions) and removes 'Virtual' atoms if
    no virtual atoms are defined.

    Parameters
    ----------
    bonded : dict
        A dictionary of bonded interactions for a single molecule,
        typically generated by `_parse_bonded`.

    Returns
    -------
    dict
        A filtered dictionary where empty interaction types and terms
        have been removed.
    """
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
    """Remove empty and 'COU' (Coulomb) interactions from a nonbonded dictionary.

    This function processes a nonbonded interactions dictionary, filtering
    out any entries that are empty or represent Coulomb ('COU') interactions.

    Parameters
    ----------
    nonbonded : dict
        A dictionary of nonbonded interactions, typically `self.nonbonded`.

    Returns
    -------
    dict
        A new dictionary identical to the input `nonbonded` dictionary
        but with empty interaction types and all 'COU' interactions removed.
    """    cleaned_dict = {}
    for pair, int_dict in nonbonded.items():
        pair_dict = {k: v for k, v, in int_dict.items() if v}
        if 'COU' in pair_dict:
            del pair_dict['COU']
        if pair_dict:
            cleaned_dict[pair] = pair_dict
    return cleaned_dict

def _remove_netf_torq_atname(all_atoms):
    """Remove 'NETF' and 'TORQ' atom names from a list of all atoms.

    This function filters a dictionary of atom definitions, returning
    only the atom names that are not identified as 'NETF' (Net Force)
    or 'TORQ' (Torque) atoms.

    Parameters
    ----------
    all_atoms : dict
        A dictionary where keys are atom numbers and values are tuples
        `('AtomType', 'AtomName')`, typically sourced from `ATO['All']`.

    Returns
    -------
    list of str
        A list containing the atom names after filtering out 'NETF' and 'TORQ'.
    """
    all_atoms_list = []
    all_atoms_list = [v[1] for k, v in all_atoms.items() if v[1] != "NETF" and v[1] != "TORQ"]
    return all_atoms_list

def _remove_netf_torq_atnum(all_atoms):
    """Remove 'NETF' and 'TORQ' atom numbers from a list of all atoms.

    This function filters a dictionary of atom definitions, returning
    only the atom numbers that are not identified as 'NETF' (Net Force)
    or 'TORQ' (Torque) atoms.

    Parameters
    ----------
    all_atoms : dict
        A dictionary where keys are atom numbers and values are tuples
        `('AtomType', 'AtomName')`, typically sourced from `ATO['All']`.

    Returns
    -------
    list of int
        A list containing the atom numbers after filtering out 'NETF' and 'TORQ'.
    """
    removed_netf_torq = [int(k) for k, v in all_atoms.items() if v[1] != 'NETF' and v[1] !='TORQ']

    return removed_netf_torq
