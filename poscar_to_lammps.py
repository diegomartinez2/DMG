#!/usr/bin/env python
import numpy as np
from collections import OrderedDict

# Dictionary of atomic masses (in g/mol) for common elements
ATOMIC_MASSES = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.941, 'Be': 9.0122, 'B': 10.811,
    'C': 12.01, 'N': 14.0067, 'O': 15.9994, 'F': 18.9984, 'Ne': 20.1797,
    'Na': 22.9897, 'Mg': 24.305, 'Al': 26.9815, 'Si': 28.0855, 'P': 30.9738,
    'S': 32.06, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
    'Sc': 44.9559, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938,
    'Fe': 55.845, 'Co': 58.9332, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.64, 'As': 74.9216, 'Se': 78.96, 'Br': 79.904,
    'Kr': 83.798
}

def read_poscar(file_path):
    """Read a POSCAR file and extract lattice vectors, atom types, counts, and coordinates."""
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Remove empty lines and strip whitespace
    lines = [line.strip() for line in lines if line.strip()]

    # Extract title (line 1)
    title = lines[0]

    # Extract scaling factor (line 2)
    scaling_factor = float(lines[1])

    # Extract lattice vectors (lines 3-5)
    lattice_vectors = np.array([
        [float(x) for x in lines[i].split()] for i in range(2, 5)
    ]) * scaling_factor

    # Extract atom types and counts (VASP 5 format: line 6 for symbols, line 7 for counts)
    atom_symbols = lines[5].split()
    atom_counts = [int(x) for x in lines[6].split()]
    total_atoms = sum(atom_counts)

    # Determine coordinate type (line 8: 'Direct' or 'Cartesian')
    coord_type = lines[7].strip().lower()
    if coord_type not in ['direct', 'cartesian']:
        raise ValueError("Coordinate type must be 'Direct' or 'Cartesian'")

    # Read coordinates (lines 9 to 9+total_atoms)
    coords = np.array([
        [float(x) for x in lines[i].split()[:3]] for i in range(8, 8 + total_atoms)
    ])

    # Assign atom types based on order
    atom_types = []
    for symbol, count in zip(atom_symbols, atom_counts):
        atom_types.extend([symbol] * count)

    return title, lattice_vectors, atom_symbols, atom_counts, coord_type, coords, atom_types

def fractional_to_cartesian(coords, lattice_vectors):
    """Convert fractional coordinates to Cartesian using lattice vectors."""
    return np.dot(coords, lattice_vectors)

def compute_lammps_box(lattice_vectors):
    """Compute LAMMPS box parameters (xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz)."""
    a, b, c = lattice_vectors

    # LAMMPS box convention:
    # a = (ax, 0, 0)
    # b = (bx, by, 0)
    # c = (cx, cy, cz)
    xlo, xhi = 0.0, a[0]
    ylo, yhi = 0.0, b[1]
    zlo, zhi = 0.0, c[2]
    xy, xz, yz = b[0], c[0], c[1]

    return xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz

def write_lammps_data(file_path, title, lattice_vectors, atom_symbols, atom_counts, coord_type, coords, atom_types):
    """Write a LAMMPS data file."""
    total_atoms = sum(atom_counts)
    num_atom_types = len(set(atom_symbols))

    # Convert coordinates to Cartesian if needed
    if coord_type.lower() == 'direct':
        coords = fractional_to_cartesian(coords, lattice_vectors)
    # If Cartesian, coordinates are already in Cartesian (scaled by scaling factor in read_poscar)

    # Compute LAMMPS box parameters
    xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz = compute_lammps_box(lattice_vectors)

    # Map atom symbols to LAMMPS atom types (1-based indexing)
    unique_symbols = list(OrderedDict.fromkeys(atom_symbols))  # Preserve order
    symbol_to_type = {symbol: i + 1 for i, symbol in enumerate(unique_symbols)}

    with open(file_path, 'w') as f:
        # Write header
        f.write(f"# LAMMPS data file for {title}\n\n")
        f.write(f"{total_atoms} atoms\n")
        f.write(f"{num_atom_types} atom types\n\n")

        # Write box dimensions
        f.write(f"{xlo:.10f} {xhi:.10f} xlo xhi\n")
        f.write(f"{ylo:.10f} {yhi:.10f} ylo yhi\n")
        f.write(f"{zlo:.10f} {zhi:.10f} zlo zhi\n")
        f.write(f"{xy:.10f} {xz:.10f} {yz:.10f} xy xz yz\n\n")

        # Write masses
        f.write("Masses\n\n")
        for i, symbol in enumerate(unique_symbols, 1):
            mass = ATOMIC_MASSES.get(symbol, 1.0)  # Default mass 1.0 if not found
            f.write(f"{i} {mass:.4f}  # {symbol}\n")
        f.write("\n")

        # Write atoms
        f.write("Atoms\n\n")
        for i, (symbol, coord) in enumerate(zip(atom_types, coords), 1):
            atom_type = symbol_to_type[symbol]
            f.write(f"{i} {atom_type} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}  # {symbol}\n")

def convert_poscar_to_lammps(poscar_file, lammps_file):
    """Main function to convert POSCAR to LAMMPS data file."""
    try:
        # Read POSCAR
        title, lattice_vectors, atom_symbols, atom_counts, coord_type, coords, atom_types = read_poscar(poscar_file)

        # Write LAMMPS data file
        write_lammps_data(lammps_file, title, lattice_vectors, atom_symbols, atom_counts, coord_type, coords, atom_types)

        print(f"Successfully converted {poscar_file} to {lammps_file}")
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    # Example usage
    poscar_file = "POSCAR"
    lammps_file = "data.lammps"
    convert_poscar_to_lammps(poscar_file, lammps_file)
