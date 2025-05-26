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
}

def poscar_to_lammps(input_file, output_file):
    """
    Convert VASP POSCAR to LAMMPS data file in 'atomic' style using ASE.
    Compatible with older ASE versions without style='atomic'.
    Includes atomic masses and supports non-orthorhombic cells.
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
    
    # Compute LAMMPS box parameters (for triclinic cells)
    xlo, ylo, zlo = 0.0, 0.0, 0.0  # Assume origin at (0,0,0)
    xhi = cell[0][0]
    yhi = cell[1][1]
    zhi = cell[2][2]
    xy = cell[1][0]
    xz = cell[2][0]
    yz = cell[2][1]
    
    # Write LAMMPS data file
    with open(output_file, 'w') as f:
        f.write('# LAMMPS data file generated from POSCAR\n\n')
        f.write(f'{natoms} atoms\n')
        f.write(f'{ntypes} atom types\n\n')
        
        # Write cell boundaries (including tilt factors for non-orthorhombic cells)
        f.write(f'{xlo:.6f} {xhi:.6f} xlo xhi\n')
        f.write(f'{ylo:.6f} {yhi:.6f} ylo yhi\n')
        f.write(f'{zlo:.6f} {zhi:.6f} zlo zhi\n')
        if any(abs(v) > 1e-6 for v in [xy, xz, yz]):  # Only write tilt factors if non-zero
            f.write(f'{xy:.6f} {xz:.6f} {yz:.6f} xy xz yz\n')
        f.write('\n')
        
        # Write Masses section
        f.write('Masses\n\n')
        for symbol, type_id in sorted(symbol_to_type.items(), key=lambda x: x[1]):
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
    Supports non-orthorhombic cells and non-zero origins.
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
    xlo = ylo = zlo = xhi = yhi = zhi = xy = xz = yz = 0.0
    for line in lines:
        if 'xlo xhi' in line:
            xlo, xhi = map(float, line.split()[:2])
        if 'ylo yhi' in line:
            ylo, yhi = map(float, line.split()[:2])
        if 'zlo zhi' in line:
            zlo, zhi = map(float, line.split()[:2])
        if 'xy xz yz' in line:
            xy, xz, yz = map(float, line.split()[:3])
    
    # Construct cell matrix (triclinic support)
    cell = [
        [xhi - xlo, 0.0, 0.0],
        [xy, yhi - ylo, 0.0],
        [xz, yz, zhi - zlo]
    ]
    
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
                # Shift coordinates to account for non-zero origin
                x, y, z = float(parts[2]) - xlo, float(parts[3]) - ylo, float(parts[4]) - zlo
                positions.append([x, y, z])
    
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