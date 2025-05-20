#!/usr/local/bin/python
from ase.io import read, write
atoms = read('structure.dat', format='lammps-data', style='atomic')
write('POSCAR', atoms, format='vasp')
