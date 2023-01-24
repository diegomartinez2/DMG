import pandas as pd
import re
from glob import glob
import numpy as np

POPULATION = 2

def en_processing(pop_id=POPULATION):
    en_path = "vasp/energies"
    energies = pd.read_csv(en_path, engine='python', sep="  ",header=None)
    print(energies)
    energies = energies[3]*0.0734986176
    processed_energy_path= f"data/energies_supercell_population{pop_id}.dat"
    e=open(processed_energy_path,"w+")
    energies.to_csv(e, index=False, header=False, float_format="%16.12f", sep='\t', mode="w+")
    e.close()


def f_processing(parsed_args,atom_occurrencies,type_atoms):
    pop_id = parsed_args[0][0]
    n = parsed_args[0][1]
    all_files = glob("vasp/forces/*")
    all_files.sort(key=lambda f: int(re.sub('\D', '', f)))
    a = 0
    for file in all_files:
        a = a + 1
        print(file)
        forces = pd.read_csv(file, engine='python', sep="\s+", skiprows=1, header=None)
        forces.drop([0, 4], axis=1, inplace=True)
        forces[3] = pd.to_numeric(forces[3], downcast="float")
        forces = forces * 0.0388937935
        processed_forces_path = f"data/forces_population{pop_id}_" + str(a) + ".dat"
        f = open(processed_forces_path, "w+")
        N = int(n)*int(n)*int(n)
        for i in range(0, N):
            for j in range(0, len(atom_occurrencies)):
                for k in range(0,atom_occurrencies[j]):
                    force = forces.iloc[N * int(np.sum(atom_occurrencies[0:j])-atom_occurrencies[j]) + atom_occurrencies[j]*i + k]
                    force = force.to_frame()
                    force = force.transpose()
                    force.to_csv(f, index=False, header=False, float_format="%16.12f", sep='\t', mode="w+")
############### everything alright, maybe read the dynamical matrix file instead of the POSCAR!
