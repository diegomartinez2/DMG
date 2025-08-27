#!/usr/bin/env python
from ase.io import read, write
from ase.build import make_supercell
import numpy as np

# Lee el archivo POSCAR original (celda unidad de 8 átomos)
atoms_unitcell = read('POSCAR')

# Define la matriz de la supercelda (2x2x2)
supercell_matrix = [[2, 0, 0],
                    [0, 2, 0],
                    [0, 0, 2]]

# Construye la supercelda
atoms_supercell = make_supercell(atoms_unitcell, supercell_matrix)

# Escribe la supercelda en un archivo SPOSCAR (formato VASP)
write('SPOSCAR', atoms_supercell, format='vasp')

print(f"Supercelda generada con {len(atoms_supercell)} átomos y guardada en SPOSCAR.")
