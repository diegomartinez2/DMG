#!/usr/local/bin/python
import numpy as np
from ase.io import read, write
from ase import Atoms

def poscar_to_lammps(input_file, output_file):
    """
    Convert VASP POSCAR to LAMMPS data file in 'atomic' style using ASE.
    Compatible with older ASE versions without style='atomic'.
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
        f.write('0.0 {:.6f} zlo zhi\n'.format(cell[2][2]))
        f.write('\nAtoms # atomic\n\n')

        # Write atom data: atom-ID, atom-type, x, y, z
        for i in range(natoms):
            f.write(f'{i+1} {atom_types[i]} {positions[i][0]:.6f} {positions[i][1]:.6f} {positions[i][2]:.6f}\n')

    print(f"Converted {input_file} to {output_file} (LAMMPS atomic style)")

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
