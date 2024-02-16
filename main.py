import afmtogmx

base_off = r"C:\Users\rjweldon\Desktop\Wang_Laptop_College\Projects\AFMTools_v2_update\test\compare\base.off"
compare_off = r"C:\Users\rjweldon\Desktop\Wang_Laptop_College\Projects\AFMTools_v2_update\test\compare\compare.off"

base_off = afmtogmx.ReadOFF(off_loc=base_off)
compare_off = afmtogmx.ReadOFF(off_loc=compare_off)


name_translation = {'OQM' : 'O', 'EQM' : 'E', 'HQM' : 'H'}

print(
    afmtogmx.compare.gen_difference_string(base_off = base_off, compare_off = compare_off, name_translation = name_translation))