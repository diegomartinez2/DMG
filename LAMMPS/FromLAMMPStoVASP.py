#!/usr/local/bin/python
from ase.io import read, write
atoms = read('structure_lammps.dat', format='lammps-data', style='atomic')
write('POSCAR', atoms, format='vasp')
