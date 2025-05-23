#!/usr/local/bin/python

from ase.io import read, write
from ase import Atoms
import glob

def fromVASPtoLAMMPS(arg):
    atoms = read(arg, format='vasp')
    write(f'{arg}.lammps', atoms, format='lammps-data', style='atomic')
    pass

def fname(arg):
    # Example: Assign atom types
    atoms = Atoms('Si2O', positions=[[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    atoms.set_tags([1, 1, 2])  # Assign atom types (e.g., 1 for Si, 2 for O)
    write(f'{poscar}.lammps', atoms, format='lammps-data')
    pass

if __name__ == '__main__':
    print(atoms.get_chemical_symbols())  # Check atom types
    print(atoms.get_positions())         # Check positions
    print(atoms.get_initial_charges())   # Check charges, if applicable

    for poscar in glob.glob('POSCAR-*'):
        atoms = read(poscar, format='vasp')
        write(f'{poscar}.lammps', atoms, format='lammps-data', style='atomic')

def image(atoms):
    from ase.visualize import view
    write('image.png', atoms)
