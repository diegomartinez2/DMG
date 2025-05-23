#!/usr/local/bin/python

from ase.io import read, write
import glob
for poscar in glob.glob('POSCAR-*'):
    atoms = read(poscar, format='vasp')
    write(f'{poscar}.lammps', atoms, format='lammps-data', style='atomic')

def fromVASPtoLAMMPS(arg):
    atoms = read(arg, format='vasp')
    write(f'{arg}.lammps', atoms, format='lammps-data', style='atomic')
    pass
