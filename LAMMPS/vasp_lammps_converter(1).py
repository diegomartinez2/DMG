#!/usr/bin/env python
"""
This code uses ASE to read and write data files of the atomic positions.
It can read LAMMPS and convert to VASP (POSCAR), it requires the atom types as data extra.
It can read VASP (POSCAR) and write into LAMMPS data file in 'atomic' format
This script is made so it will work with older versions of ASE.
"""
import numpy as np
from ase.io import read, write
from ase import Atoms

# Dictionary of atomic masses (in amu) for common elements
ATOMIC_MASSES = {
    'H': 1.00794, 'He': 4.00260, 'Li': 6.941, 'Be': 9.01218, 'B': 10.811,
    'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'F': 18.9984, 'Ne': 20.1797,
    'Na': 22.9897, 'Mg': 24.305, 'Al': 26.9815, 'Si': 28.0855, 'P': 30.9738,
    'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
    'Sc': 44.9559, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938,
    'Fe': 55.845, 'Co': 58.9332, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38
    # Add more elements as needed
}

def poscar_to_lammps(input_file, output_file):
    """
    Convert VASP POSCAR to LAMMPS data file in 'atomic' style using ASE.
    Compatible with older ASE versions without style='atomic'.
    Includes atomic masses for each atom type.
    """
    # Read POSCAR file
    atoms = read(input_file, format='vasp')

    # Get cell parameters and atom information
    cell = atoms.get_cell()
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()

    # Map unique chemical symbols to atom types (1-based for LAMMPS)
    unique_symbols = sorted(set(symbols))
    symbol_to_type = {symbol: i+1 for i, symbol in enumerate(unique_symbols)}
    atom_types = [symbol_to_type[symbol] for symbol in symbols]

    # Check if all symbols have defined masses
    for symbol in unique_symbols:
        if symbol not in ATOMIC_MASSES:
            raise ValueError(f"Atomic mass for element {symbol} not defined in ATOMIC_MASSES")

    # Count number of atoms and atom types
    natoms = len(atoms)
    ntypes = len(unique_symbols)

    # Write LAMMPS data file
    with open(output_file, 'w') as f:
        f.write('# LAMMPS data file generated from POSCAR\n\n')
        f.write(f'{natoms} atoms\n')
        f.write(f'{ntypes} atom types\n\n')

        # Write cell boundaries (assuming orthorhombic cell for simplicity)
        f.write('0.0 {:.6f} xlo xhi\n'.format(cell[0][0]))
        f.write('0.0 {:.6f} ylo yhi\n'.format(cell[1][1]))
        f.write('0.0 {:.6f} zlo zhi\n\n'.format(cell[2][2]))

        # Write Masses section
        f.write('Masses\n\n')
        for symbol, type_id in symbol_to_type.items():
            f.write(f'{type_id} {ATOMIC_MASSES[symbol]:.6f} # {symbol}\n')

        # Write Atoms section
        f.write('\nAtoms # atomic\n\n')
        for i in range(natoms):
            f.write(f'{i+1} {atom_types[i]} {positions[i][0]:.6f} {positions[i][1]:.6f} {positions[i][2]:.6f}\n')

    print(f"Converted {input_file} to {output_file} (LAMMPS atomic style with masses)")

def lammps_to_poscar(input_file, output_file, atom_symbols):
    """
    Convert LAMMPS data file in 'atomic' style to VASP POSCAR using ASE.
    atom_symbols is a list mapping LAMMPS atom types (1-based) to chemical symbols.
    """
    # Read LAMMPS data file manually
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Parse header to get number of atoms and atom types
    natoms = 0
    ntypes = 0
    for line in lines:
        if 'atoms' in line and not 'atom types' in line:
            natoms = int(line.split()[0])
        if 'atom types' in line:
            ntypes = int(line.split()[0])

    # Parse cell boundaries
    cell = []
    for line in lines:
        if 'xlo xhi' in line:
            xlo, xhi = map(float, line.split()[:2])
            cell.append([xhi - xlo, 0.0, 0.0])
        if 'ylo yhi' in line:
            ylo, yhi = map(float, line.split()[:2])
            cell.append([0.0, yhi - ylo, 0.0])
        if 'zlo zhi' in line:
            zlo, zhi = map(float, line.split()[:2])
            cell.append([0.0, 0.0, zhi - zlo])

    # Parse atom data
    positions = []
    atom_types = []
    reading_atoms = False
    for line in lines:
        if 'Atoms' in line:
            reading_atoms = True
            continue
        if reading_atoms and line.strip():
            parts = line.split()
            if len(parts) >= 5:  # Expect atom-ID, atom-type, x, y, z
                atom_types.append(int(parts[1]))
                positions.append([float(parts[2]), float(parts[3]), float(parts[4])])

    # Convert atom types to chemical symbols
    if len(atom_symbols) < ntypes:
        raise ValueError(f"Provided {len(atom_symbols)} atom symbols, but {ntypes} atom types found")
    symbols = [atom_symbols[t-1] for t in atom_types]

    # Create ASE Atoms object
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=[True, True, True])

    # Write POSCAR file
    write(output_file, atoms, format='vasp', vasp5=True, direct=True)

    print(f"Converted {input_file} to {output_file} (VASP POSCAR)")

# Example usage
if __name__ == "__main__":
    # Example: Convert POSCAR to LAMMPS
    poscar_to_lammps('POSCAR', 'data.lammps')

    # Example: Convert LAMMPS back to POSCAR
    # Need to specify the chemical symbols corresponding to LAMMPS atom types
    # e.g., if atom type 1 is Si and type 2 is O
    atom_symbols = ['Si', 'O']
    lammps_to_poscar('data.lammps', 'POSCAR_new', atom_symbols)
