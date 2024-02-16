"""This module contains various functions related to comparing two different force fields, denoted as
base_off and compare_off"""
def gen_difference_string(base_off, compare_off, name_translation = {}):
    """Generate title which details all nonbonded differences between compare_off and base_off

    :param base_off: ReadOFF object
    :param compare_off: ReadOFF object
    :param name_translation: dictionary, optional
    :return: String stating differences between force fields
    """
    # Compare nonbonded, exclude charges

    def clean_nonbonded(nonbonded, name_translation = {}):
        cleaned = {}
        for pair, int_dict in nonbonded.items():
            if len(int_dict) == 1 and 'COU' in int_dict:
                pass
            else:
                if pair[0] in name_translation:
                    at1 = name_translation[pair[0]]
                else:
                    at1 = pair[0]
                if pair[1] in name_translation:
                    at2 = name_translation[pair[1]]
                else:
                    at2 = pair[1]
                new_pair = (at1, at2)
                cleaned[new_pair] = {}
                for interaction in int_dict:
                    if interaction == 'COU':
                        pass
                    else:
                        cleaned[new_pair][interaction] = int_dict[interaction]
        return cleaned
    cleaned_base = clean_nonbonded(base_off.nonbonded, name_translation)
    cleaned_compare = clean_nonbonded(compare_off.nonbonded, name_translation)

    compare_string = "Base "

    # Items that were removed from 'Base' force field
    for pair, int_dict in cleaned_base.items():
        for interaction in int_dict:
            if pair not in cleaned_compare:
                if "Remove" not in compare_string:
                    compare_string += "Remove "
                compare_string += f"'{pair[0]}-{pair[1]} {interaction}' "
            else:
                if interaction not in cleaned_compare[pair]:
                    if "Remove" not in compare_string:
                        compare_string += "Remove "
                    compare_string += f"'{pair[0]}-{pair[1]} {interaction}' "

    # Items that were added to 'Base' force field
    for pair, int_dict in cleaned_compare.items():
        for interaction in int_dict:
            if pair not in cleaned_base:
                if "Add" not in compare_string:
                    compare_string += "Add "
                compare_string += f"'{pair[0]}-{pair[1]} {interaction}' "
            else:
                if interaction not in cleaned_base[pair]:
                    if "Add" not in compare_string:
                        compare_string += "Add "
                    compare_string += f"'{pair[0]}-{pair[1]} {interaction}' "

    return compare_string
