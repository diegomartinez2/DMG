#!/usr/local/bin/python

from ase.io import read, write
import glob

def fromVASPtoLAMMPS(arg):
    atoms = read(arg, format='vasp')
    write(f'{arg}.lammps', atoms, format='lammps-data', style='atomic')
    pass

if __name__ == '__main__':
    for poscar in glob.glob('POSCAR-*'):
        atoms = read(poscar, format='vasp')
        write(f'{poscar}.lammps', atoms, format='lammps-data', style='atomic')
