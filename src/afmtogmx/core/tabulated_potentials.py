"""This class contains functions mainly related to generation of tabulated potentials, both 
bonded and nonbonded
"""
import numpy as np
from copy import deepcopy
from afmtogmx.core import functions
from afmtogmx.core.functions import _remove_empty_and_cou_interactions_nonbonded


def _gen_included_atoms(incl_mol, bonded):
    """Generates list of atoms to include when writing tabulated potentials based on molnames found in list incl_mol.

    :param incl_mol: list of molnames which tabulated potentials should be written for
    :param bonded: self.bonded
    :return: list of all atoms for which pairs selected from the list should be used to produce tabulated potentials.
    Any atoms not included in the list will not have tabulated potentials written for them, even if they have
    interactions with atoms on the list.
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


def _filter_nonbonded(nonbonded, excluded_int, incl_atoms):
    """Takes in nonbonded, and removes excluded interactions, excluded molnames, and reformats to match
    format of INTERACTION_POWER for interactions including POW, PEX, DPO, SRD. All other interaction names are truncated
    at 3 characters.

    :param nonbonded: self.nonbonded
    :param excluded_int: list containing excluded interactions; designed with INTRA/INTER interactions in mind. For
    example, if there are C1-C1 intermolecular and C1-C1 intramolecular EXP interactions, only the inter C1-C1 should
    have a tabulated potential (as the C1-C1 intra will rely on a tablep file). If the intra EXP is named as
    'EXPINTRA', then excluded_int = ['EXPINTRA'] will not write tabulated potentials for C1-C1 'EXPINTRA'
    :param incl_atoms: list containing all atoms which pair interactions may be written for
    :return: dictionary with all pairs and interactions that should have tabulated potentials written for them
    """
    variable_power_interactions = ['POW', 'PEX', 'DPO', 'SRD']  # list for which name of interaction should be modified
    #  to include the ^power term
    output_nonbonded = {}  # output dictionary
    temp_nonbonded = _remove_empty_and_cou_interactions_nonbonded(nonbonded)  # preliminary filter of nonbonded

    temp_2_nonbonded = {}
    # Remove atoms not in incl_atoms, and interactions not in excluded_int; otherwise leave the same
    for pair, int_dict in temp_nonbonded.items():
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
    interact = interaction[:3]

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


def gen_bonded_tabpam(bond_dict, spacing, length, num_tables):
    """Generate bonded tabulated parameters in format {(parameters) : [table_number, x, U(x), F(x)], ... }

    :param bond_dict: Dictionary, self.bonded[molname]
    :param spacing: float, spacing between x-values in table
    :param length: float, total length of table
    :param num_tables: int, total number of tables generated already
    :return: dictionary with completed tabulated potentials for each 'QUA' bonded potential in bond_dict
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
    """Creates directory to write to.

    :param write_to_dir: string, either path of directory to write to, or emtpy string
    :return: string containing location of directory to write to
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
    """Writes an empty nonbonded tabulated potential file, 5 nm, spacing = 0.0005 nm"""
    x_values = np.arange(0, 5.0005, 0.0005)
    zeros = np.zeros(len(x_values))
    final_zeros_tabpot = np.array([x_values, zeros, zeros, zeros, zeros, zeros, zeros]).T
    to_dir_prefix = write_to + "/" + prefix
    np.savetxt(f'{to_dir_prefix}.xvg', X = final_zeros_tabpot, fmt = '% 20.8E', newline = '\n')


def _write_bonded_tabpot(tabpot, prefix, write_to):
    col1, col2, col3 = tabpot[1], tabpot[2], tabpot[3]
    final_tabpot = np.array([col1, col2, col3]).T
    to_dir_prefix = write_to + "/" + prefix
    np.savetxt(f"{to_dir_prefix}_b{int(tabpot[0])}.xvg", X = final_tabpot, fmt = '% 20.8E', newline = '\n')

def _translate_pairs(atom_pair, name_translation):
    """This will return translated pairs of atoms as a tuple"""
    at1, at2 = atom_pair[:]
    if at1 in name_translation:
        at1 = name_translation[at1]
    if at2 in name_translation:
        at2 = name_translation[at2]
    return (at1, at2)