import numpy as np
from collections import defaultdict
from afmtogmx.core import residues


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


def _normalize_charges_old(normalization, charges, bonded, tolerance):
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


def _normalize_charges(normalization, charges, bonded, residue_dict, residue_priority):
    """Begins process of normalizing charges, properly formats required objects, returns fully neutralized self.charges
    dictionary

    :param normalization: str, unused, method of normalization
    :param charges: dictionary, preliminary non-normalized charges dictionary
    :param bonded: self.bonded
    :param residue_dict: self.residues
    :param residue_priority: dictionary, residue_priority
    """
    updated_charges = _gen_empty_charge_dict(bonded)
    fixed_atoms = _gen_preliminary_fixed_atoms(bonded, charges)
    residue_priority = residues._generate_residue_priority(residue_dict, bonded, residue_priority)
    for molname, residue_tuple in residue_priority.items():
        for resname in residue_tuple:
            _neutralize_residue(updated_charges, residue_dict, resname, charges, fixed_atoms, molname, normalization)

def _gen_preliminary_fixed_atoms(bonded, charges):
    """Produces preliminary fixed_atoms dictionary, marking atoms that are supposed to have a 0.0 charge as fixed

    :param bonded: self.bonded
    :param charges: preliminary, non-normalized charges dictionary
    :return: fixed_charges, dictionary, properly formatted
    """
    fixed_atoms = _gen_empty_charge_dict(bonded)

    for molname, charges_dict in charges.items():
        for atname, charge in charges_dict.items():
            # Fix atoms that should have a zero charge, otherwise, mark as not fixed
            fixed_atoms[molname][atname] = False if charge != 0.0 else True
    return fixed_atoms

def _neutralize_residue(updated_charges, residue_dict, resname, charges, fixed_atoms, molname, normalization):
    atoms_in_residue = residue_dict['Definitions'][molname][resname]
    modifiable_atoms = set([i for i in atoms_in_residue if not fixed_atoms[molname][i]])
    hold_charges = {}
    # set current resname charge dict
    print(resname)
    for atom in set(atoms_in_residue):
        if atom in modifiable_atoms:
            hold_charges[atom] = round(charges[molname][atom], 5)
        else:
            hold_charges[atom] = charges[molname][atom]

    total_charge = 0.0
    for atom in atoms_in_residue:
        total_charge += hold_charges[atom]
    print(total_charge)



    pass

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
    atom_instances = defaultdict(int)  # count number of atoms per molecule
    abs_total_charge = 0
    for atom in all_atoms:
        atom_instances[atom] += 1
        abs_total_charge += abs(charges[atom])

    atom_instances = dict(sorted(atom_instances.items(), key=lambda item: item[1], reverse=True))  # sort number of
    # atoms per molecule as from highest number to lowest number of atoms
    if normalization == "M-POPULOUS":
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

    for atom, charge in charges.items():
        charges[atom] = round(charge, 5)

    return charges  # return neutralized dictionary
