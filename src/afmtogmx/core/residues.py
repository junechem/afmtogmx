def _check_residue_definitions(bonded, residue_definition):
    """Make sure that all residues and atom types are present in the .off file

    :param bonded: self.bonded
    :param residue_definition: dictionary
    """
    for molname, residues in residue_definition.items():
        if molname not in bonded:
            print(f"Molname {molname} from residue_definitions not found in .off file. Check input")
            exit(1)
        coulomb_types = set(
            [v[1] for k, v in bonded[molname]['ATO']['All'].items() if v[1] != 'NETF' and v[1] != 'TORQ'])
        vdw_types = set(
            [v[0] for k, v in bonded[molname]['ATO']['All'].items() if v[0] != 'NETF' and v[0] != 'TORQ'])
        all_atomtypes = coulomb_types.union(vdw_types)

        for residue, attype_list in residues.items():
            residue_attypes = set(attype_list)
            if not residue_attypes.issubset(all_atomtypes):
                not_in_off = residue_attypes.difference(all_atomtypes)
                print(f"Atoms {not_in_off} defined in residue_definition are not found in molname {molname}. Check your"
                      f" input")
                exit(1)


def _check_residue_atnums(bonded, residue_atnums):
    """This is a check to make sure that atom numbers specified in residue_atnums are actually in the .off file

    :param bonded: self.bonded
    :param residue_atnums: residue_atnums dictionary
    """
    for molname, residues in residue_atnums.items():
        if molname not in bonded:
            print(f"Molname {molname} from residue_definitions not found in .off file. Check input")
            exit(1)
        molname_atnums = set([int(i) for i in bonded[molname]['ATO']['All'].keys()])

        all_atnums = set()
        for residue, atnum_list in residues.items():
            for atnums in atnum_list:
                all_atnums = all_atnums.union(set(atnums))
        if not all_atnums.issubset(molname_atnums):
            not_in_off = all_atnums.difference(molname_atnums)
            print(f"Atom numbers {not_in_off} defined in residue_atnums are not found in molname {molname}. Check your "
                  f"input files")
            exit(1)


def _check_molname_resname(bonded, residues, residue_priority):
    """Checks to see if molnames, resnames defined in residue_priority are present in self.residues

    :param bonded: self.bonded
    :param residues: self.residues
    :param residue_priority: dictionary
    """
    for molname, priorities in residue_priority.items():
        if molname not in bonded:
            print(f"\nMolname {molname} from residue_priority not found in .off file. Please check residue_priority "
                  f"input")
            exit(1)
        for resname in priorities:
            if resname not in residues['Definitions'][molname] or resname not in residues['Residues'][molname]:
                print(f"\nResidue name {resname} not found in defined off.residues dictionary. Check residue_priority "
                      f"input")
                exit(1)
