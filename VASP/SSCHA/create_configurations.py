from __future__ import print_function
from __future__ import division
import sys,os

import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble
import numpy as np

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
nat = ens.current_dyn.structure.N_atoms * np.prod(ens.current_dyn.GetSupercell())
atm = np.unique(ens.current_dyn.structure.atoms)
ss = ens.current_dyn.structure.generate_supercell(ens.current_dyn.GetSupercell())
type_dict = {x : i for i, x in enumerate(atm)}
inv_dict = {i : x for x, i in type_dict.items()}
coordenadas = ens.xats.reshape((ens.N,nat,3))
line_types = " ".join([inv_dict[x] for x in np.arange(len(type_dict))])
#line_atoms = " ".join([str(type_dict[x]) for x in ss.atoms])
line_atoms = [(type_dict[x]) for x in ss.atoms]
#for i in range(N_RANDOM):
#    print (coordenadas[i])
#    print ("---------------")
    #bash_command =
    #os.sysyem(bash_command)
#print (coordenadas[0,0,:])
#print (line_types)
#print (line_atoms)
# for i in range(N_RANDOM):
#     a = np.insert(coordenadas[i],[0],np.array(line_atoms).reshape((len(line_atoms),1)),axis=1)
#     ind = np.argsort(a[:,1])
#     coordenadas_temp = a[ind]   #no funciona??
#     #now we write the POSCAR file for the 'i' ensamble
X=coordenadas[0]
Y=line_atoms
#print (zip(Y,X))
#Z=[x for _, X in sorted(zip(Y,X))]
Z= [x for (y,x) in sorted(zip(Y,X), key=lambda pair: pair[0])]
#print (Z)

#POSCAR data

# print ("SrTiO3")
# print ("a")
# print ("15.5785322165032696  0.0000000000000000  0.0000000000000000")
# print ("0.0000000000000000  15.5785322165032696  0.0000000000000000")
# print ("0.0000000000000000  0.0000000000000000  15.5785322165032696")
# print (line_types)
# print (line_atoms.count(0), line_atoms.count(1), line_atoms.count(2))
# print ("direct")
# for i in range(len(line_atoms)):
#     print (" ".join(str(e) for e in Z[i]))
with open('readme.txt', 'w') as f:
    f.write("SrTiO3")
    f.write('\n')
    f.write("a")
    f.write('\n')
    f.write("15.5785322165032696 0.0000000000000000 0.0000000000000000")
    f.write('\n')
    f.write("0.0000000000000000 15.5785322165032696 0.0000000000000000")
    f.write('\n')
    f.write("0.0000000000000000 0.0000000000000000 15.5785322165032696")
    f.write('\n')
    f.write(line_types)
    f.write('\n')
    f.write(str(line_atoms.count(0))+" ")
    f.write(str(line_atoms.count(1))+" ")
    f.write(str(line_atoms.count(2)))
    f.write('\n')
    f.write("direct")
    f.write('\n')
    for i in range(len(line_atoms)):
        f.write(" ".join(str(e) for e in Z[i]))
        f.write('\n')
