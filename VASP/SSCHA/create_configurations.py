from __future__ import print_function
from __future__ import division
import sys,os
import shutil

import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble
import numpy as np

# Define constants
Ev_to_Ry =  0.073498644351
Angst_to_bohr = 1.8897259886
# Define input variables

NQIRR = 10                    # The number of irreducible q points in the q-point grid
T = 50                        # The temperature at which we want to perform the calculation in Kelvin
SUPERCELL = (4,4,4)           # The size of the supercell (or the q point grid)
N_RANDOM = 128                # The number of configurations that will be created
POPULATION = 0                # The population to generate

# Load the starting dynamical matrices in the grid

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
src_path_ML_FF = r"./ML_FF"
src_path_bashexe = r"./run.bash"
shutil.copy(r"./run.sh", 'population'+str(POPULATION)+'_ensemble/run.sh')
nat = ens.current_dyn.structure.N_atoms * np.prod(ens.current_dyn.GetSupercell())
atm = np.unique(ens.current_dyn.structure.atoms)
ss = ens.current_dyn.structure.generate_supercell(ens.current_dyn.GetSupercell())
type_dict = {x : i for i, x in enumerate(atm)}
inv_dict = {i : x for x, i in type_dict.items()}
coordenadas = ens.xats.reshape((ens.N,nat,3))
line_types = " ".join([inv_dict[x] for x in np.arange(len(type_dict))])
#line_atoms = " ".join([str(type_dict[x]) for x in ss.atoms])
line_atoms = [(type_dict[x]) for x in ss.atoms]

dst_path_ML_FF = 'population'+str(POPULATION)+'_ensemble/ML_FF'
shutil.copy(src_path_ML_FF, dst_path_ML_FF)

for i in range(N_RANDOM):
    #dst_path_ML_FF = 'population'+str(POPULATION)+'_ensemble/'+str(i+1)+'/ML_FF'
    dst_path_bashexe = 'population'+str(POPULATION)+'_ensemble/'+str(i+1)+'/run.bash'

    os.makedirs('population'+str(POPULATION)+'_ensemble/'+str(i+1))
    #shutil.copy(src_path_ML_FF, dst_path_ML_FF)
    shutil.copy(src_path_bashexe, dst_path_bashexe)

#we write the KPOINTS file
    with open('population'+str(POPULATION)+'_ensemble/'+str(i+1)+'/KPOINTS', 'w') as f:
        lines = ["Not only Gamma point\n", "0\n", "Gamma\n", "2 2 2\n", "0 0 0\n" ]
        # f.write("Not only Gamma point"+'\n')
        # f.write("0"+'\n')
        # f.write("Gamma"+'\n')
        # f.write("2 2 2"+'\n')
        # f.write("0 0 0"+'\n')
        f.writelines(lines)
        f.close()

#we write the POTCAR (is this necessary for mechine learning potentials?)
    #remember that the order of the atoms in the POTCAR must be the same as the POSCAR
    shutil.copy(r"./POTCAR", 'population'+str(POPULATION)+'_ensemble/'+str(i+1)+'/POTCAR')

#we write the INCAR file (Check this out!!!)
    with open('population'+str(POPULATION)+'_ensemble/'+str(i+1)+'/INCAR', 'w') as f:
        f.write("SYSTEM = SrTiO3"+'\n')
        f.write("ISTART = 0"+'\n') # this is a new job
        f.write('\n') #ionic relaxation
        f.write("IBRION = -1"+'\n') #-1 = no ionic displacements #2 = ionic relaxation (conjugate gradient algorithm)
        #f.write("NSW    = 1") #no of ionic steps
        f.write("NELM = 25"+'\n') # maximum of 25 electronic steps (SCF cycles)
        f.write("EDIFF  =  1E-09"+'\n') #(SCF energy convergence; in eV)
        #f.write("ISIF   = 3") #update positions, cell shape and volume
        f.write("PREC   =  Accurate "+'\n') #  (Precision level) Must be accurate to get good stresses values
        f.write("ENCUT  =  600"+'\n') #(Cut-off energy for plane wave basis set, in eV) should be set explicitely in addition to AT LEAST 1.3* max(ENMAX) of the PPs. USE "grep ENMAX POTCAR"
        f.write('\n') #machine learning
        f.write("ML_LMLFF  = T"+'\n')
        f.write("ML_ISTART = 2"+'\n')
        f.write("RANDOM_SEED =         688344966                0                0"+'\n')
        f.close()
#now we write the POSCAR file for the 'i' ensamble
    X=coordenadas[i]
    Y=line_atoms
#print (zip(Y,X))
#Z=[x for _, X in sorted(zip(Y,X))]
    Z= [x for (y,x) in sorted(zip(Y,X), key=lambda pair: pair[0])]
#print (Z)

#POSCAR data

    with open('population'+str(POPULATION)+'_ensemble/'+str(i+1)+'/POSCAR', 'w') as f:
        f.write("SrTiO3")
        f.write('\n')
        f.write("1.0")
        f.write('\n')
        f.write(str(SUPERCELL[0]*ens.dyn_0.alat/Angst_to_bohr)+" 0.00000"+" 0.00000")
        f.write('\n')
        f.write("0.00000 "+str(SUPERCELL[1]*ens.dyn_0.alat/Angst_to_bohr)+" 0.00000")
        f.write('\n')
        f.write("0.00000 "+"0.00000 "+str(SUPERCELL[2]*ens.dyn_0.alat/Angst_to_bohr))
        f.write('\n')
        f.write(line_types)
        f.write('\n')
        f.write(str(line_atoms.count(0))+" ")
        f.write(str(line_atoms.count(1))+" ")
        f.write(str(line_atoms.count(2)))
        f.write('\n')
        f.write("Cartesian") # f.write("direct")
        f.write('\n')
        for j in range(len(line_atoms)):
            f.write(" ".join(str(round(e/Angst_to_bohr,6)) for e in Z[j]))  #rounding to 5 or 6 decimals because VASP fails if not rounded
            f.write('\n')
        f.close()
