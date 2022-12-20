import os
#exitcode = os.system('vasp')
exitcode = os.system('mpirun -n 16 vasp')   # to run parallel calculations
