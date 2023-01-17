from __future__ import print_function
from __future__ import division
import sys,os

import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble

# Define input variables

NQIRR = 10                    # The number of irreducible q points in the q-point grid
T = 50                        # The temperature at which we want to perform the calculation in Kelvin
SUPERCELL = (4,4,4)           # The size of the supercell (or the q point grid)
N_RANDOM = 100                # The number of configurations that will be created
POPULATION = 0                # The population to generate

# Load the starting dynamical matrices in the 2x2x2 grid

if (POPULATION == 0):
    namefile='harmonic_calculation/harmonic_dyn'  # We will start from the harmonic dynamical matrices at the first step
else:
    namefile='dyn_end_population'+str(POPULATION-1)+'_'  # We will start from the output dynamical matrices at the previous step

dyn = CC.Phonons.Phonons(namefile, NQIRR)

# Apply the sum rule and symmetries

dyn.Symmetrize()

# Flip the imaginary frequencies into real ones if there are imaginary phonon frequencies

dyn.ForcePositiveDefinite()


# We can print the frequencies to show the magic:

w_s, pols = dyn.DiagonalizeSupercell()
print ("\n".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in  w_s]))

# Generate the ensemble

ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL)
ens.generate(N_RANDOM, sobol = True, sobol_scatter = 0.0)

# Save the ensemble

namefile='population'+str(POPULATION)+'_ensemble'
ens.save(namefile, POPULATION)
#ens.save_extxyz(namefile, POPULATION)
#self.structures -> get_ase_atoms()
#ens.save_enhanced_xyz(

# Prepare qe input files

for i in range(N_RANDOM):
    bash_command = 'cat head_scf.in population'+str(POPULATION)+'_ensemble/scf_population'+str(POPULATION)+'_'+str(i+1)+'.dat > population'+str(POPULATION)+'_ensemble/scf_'+str(i+1)+'.in'
    os.system(bash_command)

# Prepare VASP input files

#for i in range(N_RANDOM):
    #bash_command =
    #os.sysyem(bash_command)
