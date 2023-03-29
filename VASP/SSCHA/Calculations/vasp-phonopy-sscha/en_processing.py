#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd

def en_processing(pop_id=1):
    en_path = "vasp/energies"
    energies = pd.read_csv(en_path, engine='python', sep="  ",header=None)
    print(energies)
#    energies = energies*0.0734986176
    processed_energy_path= f"data/energies_supercell_population{pop_id[0][0]}.dat"
    e=open(processed_energy_path,"w+")
    energies.to_csv(e, index=False, header=False, float_format="%16.12f", sep='\t', mode="w+")
    e.close()
