"""This class contains functions mainly related to generation of tabulated potentials, both 
bonded and nonbonded
"""
import copy

import numpy as np
from copy import deepcopy
from afmtogmx.core import functions
from afmtogmx.core.functions import _remove_empty_and_cou_interactions_nonbonded


def _gen_included_atoms(incl_mol, bonded):
    """Generate a list of atom names to include for tabulated potential generation.

    This function determines which atom names should be considered for generating
    tabulated potentials. If a list of molecule names (`incl_mol`) is provided,
    only atoms belonging to those molecules (and not 'NETF'/'TORQ' types)
    will be included. Otherwise, all non-'NETF'/'TORQ' atom names from the
    `bonded` data will be included.

    Parameters
    ----------
    incl_mol : list of str or None
        A list of molecule names. If provided, only atoms from these molecules
        will be included. If `None` or empty, all relevant atoms are considered.
    bonded : dict
        The `self.bonded` dictionary from the `ReadOFF` class, containing
        parsed bonded interaction information, specifically the 'ATO' section.

    Returns
    -------
    list of str
        A list of unique atom names (strings) that are designated for
        tabulated potential generation.
    """

    all_atoms = list(set([atom[0] for key, value in bonded.items() for atnum, atom in
                          value['ATO']['All'].items() if atom[0] != "NETF" and atom[0] != "TORQ"]))
    if incl_mol:  # if only including certain molnames, generate list with just those atoms
        atoms_to_remove = all_atoms
        for molname in incl_mol:
            molname_atoms = [atom[0] for atnum, atom in bonded[molname]['ATO']['All'].items() if atom[0] != "NETF" and
                             atom[0] != "TORQ"]
            atoms_to_remove = list(set(molname_atoms) ^ set(atoms_to_remove))  #
        included_atoms = list(set(all_atoms) - set(atoms_to_remove))
    else:  # otherwise, just include all atoms from .off file
        included_atoms = all_atoms

    return included_atoms


def _filter_nonbonded(nonbonded, excluded_int, incl_atoms, excl_pairs):
    """Filter and reformat nonbonded interactions for potential generation.

    This function processes a nonbonded dictionary by removing interactions
    and pairs specified for exclusion, filtering based on a list of
    included atoms, and renaming certain interaction types to a consistent
    format (e.g., 'POW' becomes 'POW_6' if the power is 6).

    Parameters
    ----------
    nonbonded : dict
        A dictionary of nonbonded interactions, typically `self.nonbonded`
        from the `ReadOFF` class.
    excluded_int : list of str
        A list of interaction names (as they appear in the .off file)
        to be excluded from the returned dictionary.
    incl_atoms : list of str
        A list of atom names that are designated for inclusion. Any
        interaction pair containing an atom not in this list will be
        excluded.
    excl_pairs : list of list
        A list of atom pairs (e.g., `[['C', 'H']]`) to be excluded from
        the returned dictionary.

    Returns
    -------
    dict
        A dictionary containing the filtered and reformatted nonbonded
        interactions that are ready for tabulated potential generation.
    """
    variable_power_interactions = ['POW', 'PEX', 'DPO', 'SRD']  # list for which name of interaction should be modified
    #  to include the ^power term
    output_nonbonded = {}  # output dictionary
    temp_nonbonded = _remove_empty_and_cou_interactions_nonbonded(nonbonded)  # preliminary filter of nonbonded

    temp_2_nonbonded = {}
    # Remove atoms not in incl_atoms, and interactions not in excluded_int; otherwise leave the same
    for pair, int_dict in temp_nonbonded.items():
        if list(pair) in excl_pairs or list(pair[::-1]) in excl_pairs:
            continue
        for interact, param_list in int_dict.items():
            if interact not in excluded_int and pair[0] in incl_atoms and pair[1] in incl_atoms:
                if pair not in temp_2_nonbonded:
                    temp_2_nonbonded[pair] = dict()
                    temp_2_nonbonded[pair][interact] = deepcopy(param_list)
                else:
                    temp_2_nonbonded[pair][interact] = deepcopy(param_list)

    for pair, int_dict in temp_2_nonbonded.items():  # populating final dictionary
        output_nonbonded[pair] = dict()

        for interaction, param_list in int_dict.items():
            if interaction[:3] in variable_power_interactions:  # Rename interactions that are variable power
                for param_set in param_list:
                    power = abs(int(param_set[1]))
                    new_int = f'{interaction[:3]}_{power}'
                    if new_int not in output_nonbonded[pair]:
                        output_nonbonded[pair][new_int] = []
                        output_nonbonded[pair][new_int].append(param_set)
                    else:
                        output_nonbonded[pair][new_int].append(param_set)
            else:  # Otherwise, just take first 3 letters of interaction
                for param_set in param_list:
                    if interaction[:3] not in output_nonbonded[pair]:
                        output_nonbonded[pair][interaction[:3]] = []
                        output_nonbonded[pair][interaction[:3]].append(param_set)
                    else:
                        output_nonbonded[pair][interaction[:3]].append(param_set)

    return output_nonbonded


def _gen_nonbond_tabpam(nonbonded, spec_pairs, spacing, length, scale_C6, charges={}):
    """Generate nonbonded tabulated potential parameters.

    This function calculates the potential and force for each nonbonded
    interaction pair over a specified range and spacing. It categorizes
    interactions into attractive and repulsive components and returns them
    in a structured dictionary.

    Parameters
    ----------
    nonbonded : dict
        A filtered dictionary of nonbonded interactions.
    spec_pairs : dict
        A dictionary defining custom attractive interactions for specific pairs.
    spacing : float
        The spacing (in nm) between x-values in the generated tables.
    length : float
        The total length (in nm) of the generated tables.
    scale_C6 : bool
        If `True`, attractive interactions are scaled by C6 for dispersion
        corrections.
    charges : dict, optional
        A dictionary of charges. Currently unused but preserved for future
        functionality.

    Returns
    -------
    dict
        A dictionary where keys are atom pairs and values are lists of NumPy
        arrays: `[x_values, COU_pot, COU_force, ATT_pot, ATT_force, REP_pot, REP_force]`.
    """
    default_attractive = ['POW_6', 'DPO_6', 'SRD_6']
    true_x_values = np.arange(0, length + spacing, spacing)
    calc_x_values = deepcopy(true_x_values)  # prohibits calc_x_values from pointer-like behavior
    calc_x_values[0] = calc_x_values[1]  # prevents problems with infinity at x = 0 for certain interactions
    nonbond_tabpot = {key: {} for key, value in nonbonded.items()}

    if spec_pairs and scale_C6:
        print("Using custom attractive functions outside of POW_6, DPO_6, PEX_6, and SRD_6 is not supported when also"
              " using scale_C6. Please rerun without scale_C6 if using custom attractive potentials")
        exit(1)

    for pair, int_dict in nonbonded.items():
        num_attr = 0
        COU_pot = np.zeros(len(true_x_values))
        COU_force = np.zeros(len(true_x_values))
        ATT_pot = np.zeros(len(true_x_values))
        ATT_force = np.zeros(len(true_x_values))
        REP_pot = np.zeros(len(true_x_values))
        REP_force = np.zeros(len(true_x_values))
        if pair not in spec_pairs:
            for interaction in int_dict:
                if interaction not in default_attractive and interaction != 'THC' and interaction != 'BUC':
                    for param_set in int_dict[interaction]:
                        cur_pot, cur_force = gen_nonbond_table(x=calc_x_values, interaction=interaction,
                                                               params=param_set, ATT=False, scale_C6=False)
                        REP_pot += cur_pot
                        REP_force += cur_force
                elif interaction == 'THC':
                    print("UNABLE TO HANDLE THC FUNCTIONS AT THIS TIME: EDIT CODE YOURSELF OR ASK RAY TO DO IT")
                elif interaction == 'BUC':
                    for param_set in int_dict[interaction]:
                        P1, P2, P3 = param_set

                        cur_pot, cur_force = gen_nonbond_table(x=calc_x_values, interaction = 'POW', params = [P2, -6],
                                                               ATT=True, scale_C6 = scale_C6)

                        ATT_pot += cur_pot
                        ATT_force += cur_force

                        cur_pot, cur_force = gen_nonbond_table(x=calc_x_values, interaction = 'EXP', params = [P1, P3],
                                                               ATT=False, scale_C6 = False)
                        REP_pot += cur_pot
                        REP_force += cur_force
                else:
                    for param_set in int_dict[interaction]:
                        num_attr += 1
                        if num_attr > 1 and scale_C6:
                            print(
                                "Cannot have greater than 1 attractive interactions and properly scale C6 for dispersion "
                                "corrections. Set 'scale_C6' to false and restart to continue generating nonbonded")
                            exit(1)
                        cur_pot, cur_force = gen_nonbond_table(x=calc_x_values, interaction=interaction,
                                                               params=param_set, ATT=True,
                                                               scale_C6=scale_C6)

                        ATT_pot += cur_pot
                        ATT_force += cur_force
            nonbond_tabpot[pair] = [true_x_values, COU_pot, COU_force, ATT_pot, ATT_force, REP_pot, REP_force]

        else:
            print(f"Custom attractive potential for pair {pair}: {spec_pairs[pair]}")
            custom_attractive = spec_pairs[pair]
            print(custom_attractive)
            for interaction in int_dict:
                print(interaction)
                if interaction not in custom_attractive and interaction != 'THC':
                    for param_set in int_dict[interaction]:
                        cur_pot, cur_force = gen_nonbond_table(x=calc_x_values, interaction=interaction,
                                                               params=param_set, ATT=False, scale_C6=False)
                        REP_pot += cur_pot
                        REP_force += cur_force
                elif interaction == 'THC':
                    print("UNABLE TO HANDLE THC FUNCTIONS AT THIS TIME: EDIT CODE YOURSELF OR ASK RAY TO DO IT")
                else:
                    for param_set in int_dict[interaction]:
                        print("PARAM SET", param_set)
                        cur_pot, cur_force = gen_nonbond_table(x=calc_x_values, interaction=interaction,
                                                               params=param_set, ATT=True, scale_C6=scale_C6)
                        ATT_pot += cur_pot
                        ATT_force += cur_force
            nonbond_tabpot[pair] = [true_x_values, COU_pot, COU_force, ATT_pot, ATT_force, REP_pot, REP_force]

    return nonbond_tabpot


def gen_nonbond_table(x, interaction, params, ATT, scale_C6):
    """Calculate potential and force for a single nonbonded interaction type.

    This function serves as a dispatcher, converting units and calling the
    appropriate potential function (e.g., `exp`, `srd`, `powe`) based on
    the interaction type.

    Parameters
    ----------
    x : numpy.ndarray
        An array of distance values at which to calculate the potential/force.
    interaction : str
        The type of interaction (e.g., 'EXP', 'POW_6'). The first three
        letters are used to select the potential function.
    params : list or tuple
        A list of parameters for the potential function.
    ATT : bool
        If `True`, indicates the interaction is attractive.
    scale_C6 : bool
        If `True` and `ATT` is `True`, the C6 parameter is scaled for
        dispersion corrections.

    Returns
    -------
    tuple of (numpy.ndarray, numpy.ndarray)
        A tuple containing the calculated potential and force arrays.
        Returns (1, 2) for unhandled interaction types.
    """
    interact = interaction[:3]
    params = copy.deepcopy(params)

    if interact == 'EXP':
        params[0], params[1] = 4.184 * params[0], 10 * params[1]
        return functions.exp(params, x)
    if interact == 'STR':
        params[0], params[1], params[2] = 4.184 * params[0] * (0.1 ** (abs(params[1]))), params[1], 0.1 * params[2]
        return functions.shtr(params, x)
    if interact == 'POW':
        params[0], params[1] = 4.184 * params[0] * (0.1 ** (abs(params[1]))), params[1]
        if ATT and scale_C6:
            params[0] = -1
        return functions.powe(params, x)
    if interact == 'SRD':
        params[0], params[1], params[2] = 4.184 * params[0] * (0.1 ** abs(params[1])), abs(params[1]), 0.1 * params[2]
        if ATT and scale_C6:
            params[0] = -1
        return functions.srd(params, x)
    if interact == 'PEX':
        print("CANNOT HANDLE PEX AT THE MOMENT: ADD FXN TO tabulated_potentials.py OR ASK RAY TO DO IT")
        return 1, 2
    if interact == 'DPO':
        print("CANNOT HANDLE DPO AT THE MOMENT: ADD FXN TO tabulated_potentials.py OR ASK RAY TO DO IT")
        return 1, 2
    else:
        print(f"Interaction {interact} not yet added to script. Add fxn to tabulated_potentials.py or ask Ray to do it")


def _scale_for_FE(sc_sigma, nonbonded_string, tabpot):
    """Scale nonbonded tabulated potentials for free energy (FE) calculations.

    This function adjusts the repulsive component of tabulated potentials
    to be compatible with free energy calculations that use the `sc-sigma`
    value in a GROMACS `grompp.mdp` file. It scales the C12 parameter based
    on the provided `sc_sigma` value.

    Parameters
    ----------
    sc_sigma : float
        The `sc-sigma` value that will be used in the `grompp.mdp` file.
    nonbonded_string : str
        The string content of the `[ nonbond_params ]` section, used to
        extract C6 and C12 values for each pair.
    tabpot : dict
        The dictionary of nonbonded tabulated potentials, as generated by
        `_gen_nonbond_tabpam`.

    Returns
    -------
    dict
        The `tabpot` dictionary with its repulsive potential and force
        components (`tabpot[pair][5]` and `tabpot[pair][6]`) scaled for
        the free energy calculation.
    """

    for line in nonbonded_string.split('\n'):
        entries = line.split()
        if len(entries) == 0:
            continue
        C6, C12 = float(entries[-2]), float(entries[-1])
        if C6 != 0 and C12 != 0:
            atom_pair = tuple(entries[0:2])
            new_c12 = C6 * (sc_sigma**6)
            scale_c12 = 1/new_c12
            tabpot[atom_pair][5] = tabpot[atom_pair][5]*scale_c12
            tabpot[atom_pair][6] = tabpot[atom_pair][6]*scale_c12
        else:
            pass
    return tabpot




def gen_bonded_tabpam(bond_dict, spacing, length, num_tables):
    """Generate bonded tabulated potential parameters for 'QUA' and 'QBB' interactions.

    This function calculates the potential and force for quartic 'QUA' bonds
    and 'QBB' bond-angle interactions over a specified range and spacing.
    It returns a dictionary containing the tabulated data.

    Parameters
    ----------
    bond_dict : dict
        A dictionary of bonded interactions for a single molecule, typically
        `self.bonded[molname]`.
    spacing : float
        The spacing (in nm) between x-values in the generated tables.
    length : float
        The total length (in nm) of the generated tables.
    num_tables : int
        The starting number to use for table indexing. This allows for
        unique table numbers across multiple molecules.

    Returns
    -------
    tuple of (dict, int)
        - `bonded_tabpam_dict`: A dictionary where keys are parameter tuples and
          values are lists `[table_number, x_values, U(x), F(x)]`.
        - `num_tables`: The updated total number of tables generated.
    """
    if not bond_dict['BON']['QUA'] and not bond_dict['BD3']['QBB']:  # If 'QUA' in bond_dict['BON'] not populated, skip
        return {}, num_tables
    else:
        true_x_values = np.arange(0, length + spacing, spacing)
        calc_x_values = deepcopy(true_x_values)  # prohibits calc_x_values from pointer-like behavior
        calc_x_values[0] = calc_x_values[1]  # prevents problems with infinity at x = 0 for certain interactions
        bonded_tabpam_dict = dict()  # dictionary to hold completed tabulated paramters
        for params, atom_list in bond_dict['BON']['QUA'].items():
            # Modify parameters to have proper units
            P1, P2, P3, P4 = 1E-1 * params[0], 4.184E2 * params[1], 4.184E3 * params[2], 4.184E4 * params[3]
            potential, force = functions.quarbond(param_list=[P1, P2, P3, P4], r=calc_x_values)  # get potential, force
            bonded_tabpam_dict[params] = [num_tables, true_x_values, potential, force]  # add to dictionary
            num_tables += 1  # keep track of the different table numbers
        for params, atom_list in bond_dict['BD3']['QBB'].items():  # Accounting for QBB terms
            P1, P2, P3, P4, P5 = 1E-1 * params[0], 4.184 * 1E2 * params[1], 4.184E2 * params[2], 4.184E3 * params[
                3], 4.184E4 * params[4]
            potential, force = functions.quarbond(param_list=[P1, P3, P4, P5], r=calc_x_values)
            bonded_tabpam_dict[params] = [num_tables, true_x_values, potential, force]
            num_tables += 1

    return bonded_tabpam_dict, num_tables  # return


def _to_dir(write_to_dir=""):
    """Ensure the existence of and return the path to an output directory.

    If a directory path is provided, this function checks if it exists,
    creating it if necessary. If no path is provided, it creates a
    default `tabpot` subdirectory in the current working directory.

    Parameters
    ----------
    write_to_dir : str, optional
        The path to the target directory. If empty, a default directory
        named `tabpot` will be created in the current working directory.

    Returns
    -------
    str
        The validated path to the output directory.
    """
    from re import sub
    from os.path import exists
    from os import getcwd, mkdir
    if write_to_dir:
        if not exists(write_to_dir):
            mkdir(write_to_dir)
        else:
            pass
        return write_to_dir
    else:
        CWD = sub(r"\\", "/", getcwd())
        write_to = f'{CWD}/tabpot/'
        if not exists(write_to):
            mkdir(write_to)
        return write_to

def _write_nonbonded_pair_tabpot(atom_pair, tabpot, name_translation, write_to, prefix):
    """Write the tabulated potential for a single nonbonded atom pair to an .xvg file.

    Parameters
    ----------
    atom_pair : tuple of str
        The pair of atom names, e.g., `('C', 'H')`.
    tabpot : list of numpy.ndarray
        A list of tabulated potential data arrays for the pair.
    name_translation : dict
        A dictionary for translating atom names for use in the filename.
    write_to : str
        The path to the output directory.
    prefix : str
        The prefix to use for the output filename. The final filename will be
        `<write_to>/<prefix>_<At1>_<At2>.xvg`.
    """
    col1, col2, col3, col4, col5, col6, col7 = tabpot #[np.reshape(tabpot[i], (x_vals, 1)) for i in range(0, 7)]
    final_tabpot = np.array([col1, col2, col3, col4, col5, col6, col7]).T

    At1 = atom_pair[0]
    At2 = atom_pair[1]
    if atom_pair[0] in name_translation:
        At1 = name_translation[atom_pair[0]]
    if atom_pair[1] in name_translation:
        At2 = name_translation[atom_pair[1]]
    to_dir_prefix = write_to + "/" +  prefix

    np.savetxt(f"{to_dir_prefix}_{At1}_{At2}.xvg", X = final_tabpot, fmt = '% 20.8E', newline = '\n')

def _write_blank_nonbonded(prefix, write_to):
    """Write an empty nonbonded tabulated potential file required by GROMACS.

    GROMACS requires a default, blank table file. This function generates
    a file of all zeros out to 5.0 nm with a spacing of 0.0005 nm.

    Parameters
    ----------
    prefix : str
        The prefix to use for the output filename. The final filename will be
        `<write_to>/<prefix>.xvg`.
    write_to : str
        The path to the output directory.
    """
    x_values = np.arange(0, 5.0005, 0.0005)
    zeros = np.zeros(len(x_values))
    final_zeros_tabpot = np.array([x_values, zeros, zeros, zeros, zeros, zeros, zeros]).T
    to_dir_prefix = write_to + "/" + prefix
    np.savetxt(f'{to_dir_prefix}.xvg', X = final_zeros_tabpot, fmt = '% 20.8E', newline = '\n')


def _write_bonded_tabpot(tabpot, prefix, write_to):
    """Write the tabulated potential for a single bonded interaction to an .xvg file.

    Parameters
    ----------
    tabpot : list
        A list of tabulated potential data arrays for the bonded interaction,
        in the format `[table_number, x_values, U(x), F(x)]`.
    prefix : str
        The prefix to use for the output filename.
    write_to : str
        The path to the output directory. The final filename will be
        `<write_to>/<prefix>_b<table_number>.xvg`.
    """
    col1, col2, col3 = tabpot[1], tabpot[2], tabpot[3]
    final_tabpot = np.array([col1, col2, col3]).T
    to_dir_prefix = write_to + "/" + prefix
    np.savetxt(f"{to_dir_prefix}_b{int(tabpot[0])}.xvg", X = final_tabpot, fmt = '% 20.8E', newline = '\n')

def _translate_pairs(atom_pair, name_translation):
    """Translate atom names within a pair using a name translation dictionary.

    Parameters
    ----------
    atom_pair : tuple of str
        The pair of atom names to be translated, e.g., `('C_off', 'H_off')`.
    name_translation : dict
        A dictionary mapping original atom names to their desired new names,
        e.g., `{'C_off': 'C_top'}`.

    Returns
    -------
    tuple of str
        A new tuple containing the translated atom names. If an atom name
        is not found in the `name_translation` dictionary, the original
        name is used.
    """
    at1, at2 = atom_pair[:]
    if at1 in name_translation:
        at1 = name_translation[at1]
    if at2 in name_translation:
        at2 = name_translation[at2]
    return (at1, at2)