"""
import os
#exitcode = os.system('vasp')
exitcode = os.system('mpirun -n 16 vasp')   # to run parallel calculations
"""



"""
    #!/usr/bin/env bash
    #SBATCH  LINES …....................

	export NCORES=$(($SLURM_JOB_NUM_NODES*16))

        cp inputfiles $TMPDIR
        cd $TMPDIR

        module purge # unloads all modules, so to have a fresh start
                               # and no competing modules
        module load iccifort/2019.1.144-GCC-8.2.0-2.31.1
        module load impi/2018.4.274
        module load Python/3.7.2
        module load ASE/3.18.0-Python-3.7.2
        module load VASP/5.4.4.16052018

	time python test_pd111.py $SLURM_SUBMIT_DIR > @SLURM_SUBMIT_DIR/result.txt

        exit_code=$?

        echo ""
        echo "Executable 'vasp' finished with exit code $exit_code"

       # Copy the files back after run:

       echo " "
       echo "Remove unneeded files"
       rm -vf $TMPDIR/{Unneeded file, comma-sperated}

       echo " "
       echo "copy needed files back"
       mv $TMPDIR/* $SLURM_SUBMIT_DIR
"""

"""
.bashrc
        export VASP_SCRIPT="$HOME/VASP/run_vasp.py"
        export VASP_PP_PATH="$HOME/VASP/"
"""

"""
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
from ase.build import fcc111

# creates a Pt(111) surface from the ase.build model with standard lattice constant
atoms=fcc111('Pt', size=(2,2,4))
calc = vasp_calculator.Vasp()      # Define whatever parameters you need
calc.initialize(atoms)         #Initialize the VASP calculator class with the atom-object
calc.write_input(atoms)    #Call the write_input definition of the VASP calculator class
"""

from ase.calculators.vasp import Vasp
from ase.optimize import QuasiNewton
from ase.constraints import FixAtoms
from ase.build import fcc111

#Setup Pd
a0 = 3.996
k0 = 6
size = [1,1,4]

syst = fcc111(symbol='Pd',size=size,a=a0,vacuum=8.0)

c = FixAtoms(indices=[atom.index for atom in syst if atom.tag<= 3])
syst.set_constraint(c)

calc = Vasp(	prec='Normal',xc='PBE',
		encut = 450.0, #PW cut-off
		ispin=1, #No spin-polarization
		ibrion=-1, #No VASP relaxation.
		nsw = 0, #Max. no of relaxation steps.
		kpts = [k0,k0,1],
		sigma=0.05,
		potim=0.5,
		isif=0,
		npar=4 #Band parallezation use kpar for k-points
		)

syst.set_calculator(calc) #Connect atoms and calculator

#Now run the relaxation with QuasiNewton
dyn = QuasiNewton(syst,trajectory='pd-relax.traj')
dyn.run(fmax=0.05)

print ("Found energy after relaxing", syst.get_potential_energy() )
