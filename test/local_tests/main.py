import afmtogmx
import os

base_off = os.path.dirname(afmtogmx.__file__) + r"/../../test/compare/base.off"
compare_off = os.path.dirname(afmtogmx.__file__) + r"/../../test/compare/compare.off"


base_off = afmtogmx.ReadOFF(off_loc=base_off)
compare_off = afmtogmx.ReadOFF(off_loc=compare_off)


name_translation = {'OQM' : 'O', 'EQM' : 'E', 'HQM' : 'H'}

print(
    afmtogmx.compare.gen_difference_string(base_off = base_off, compare_off = compare_off, name_translation = name_translation))