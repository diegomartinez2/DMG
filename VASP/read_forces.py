# This is an example about how to read forces in an OUTCAR file
#
# The forces acting on each ions are provided in an OUTCAR file.
# you can read the table using the Outcar.read_table_pattern() method
# here are the working regex in order to read the forces.
# At the end you get a multiple list with forces as float.
# With last_one_only=False, the first index of the list is the ionic step,
# the second is the atom index and the last is 0, 1 or 2 for f_x, f_y and f_z
# respectively.
#
# The script read the following table in the OUTCAR file:
#
# POSITION                                       TOTAL-FORCE (eV/Angst)                                                                                                         
# -----------------------------------------------------------------------------------
#     52.57643     30.37341      1.83720         0.007660      0.002277     -0.000107
#     40.97214     49.74181      1.83706         0.002027      0.005968     -0.002330
#     18.39530     49.36766      1.83714        -0.005405      0.005473     -0.000319
# [...]
#      7.42357     29.62659      1.83720        -0.007660     -0.002277     -0.000107
#     19.02786     10.25819      1.83706        -0.002027     -0.005968     -0.002330
#     41.60470     10.63234      1.83714         0.005405     -0.005473     -0.000319
# -----------------------------------------------------------------------------------
#

from pymatgen.io.vasp.outputs import Outcar

outcar = Outcar("OUTCAR")

forces = outcar.read_table_pattern(
    header_pattern=r"\sPOSITION\s+TOTAL-FORCE \(eV/Angst\)\n\s-+",
    row_pattern=r"\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+[+-]?\d+\.\d+\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
    footer_pattern=r"\s--+",
    postprocess=lambda x: float(x),
    last_one_only=False
)