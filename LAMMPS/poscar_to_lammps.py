import numpy as np

# Dictionary of atomic masses (amu) for common elements
ATOMIC_MASSES = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.941, 'Be': 9.0122, 'B': 10.811,
    'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'F': 18.9984, 'Ne': 20.1797,
    'Na': 22.9898, 'Mg': 24.305, 'Al': 26.9815, 'Si': 28.0855, 'P': 30.9738,
    'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
    'Sc': 44.9559, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938,
    'Fe': 55.845, 'Co': 58.9332, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.64, 'As': 74.9216, 'Se': 78.96, 'Br': 79.904,
    'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.9059, 'Zr': 91.224,
    'Nb': 92.9064, 'Mo': 95.96, 'Tc': 98.9062, 'Ru': 101.07, 'Rh': 102.9055,
    'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71,
    'Sb': 121.76, 'Te': 127.6, 'I': 126.9045, 'Xe': 131.293, 'Cs': 132.9055,
    'Ba': 137.327, 'La': 138.9055, 'Ce': 140.116, 'Pr': 140.9077, 'Nd': 144.242,
    'Pm': 144.9127, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.9254,
    'Dy': 162.5, 'Ho': 164.9303, 'Er': 167.259, 'Tm': 168.9342, 'Yb': 173.054,
    'Lu': 174.9668, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207,
    'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084, 'Au': 196.9666, 'Hg': 200.59,
    'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804, 'Po': 208.9824, 'At': 209.9871,
    'Rn': 222.0176
}

def read_poscar(filename):
    """Read a VASP POSCAR file and extract lattice, atom types, counts, and coordinates."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Line 1: Comment
    comment = lines[0].strip()
    
    # Line 2: Scaling factor
    scale = float(lines[1].strip())
    
    # Lines 3-5: Lattice vectors
    lattice = np.array([list(map(float, lines[i].split())) for i in range(2, 5)]) * scale
    
    # Line 6: Element types
    elements = lines[5].split()
    
    # Line 7: Number of atoms per type
    atom_counts = list(map(int, lines[6].split()))
    
    # Line 8: Coordinate type (Direct or Cartesian)
    coord_type = lines[7].strip().lower()
    if coord_type.startswith('s'):  # Selective dynamics
        coord_type = lines[8].strip().lower()
        start_line = 9
    else:
        start_line = 8
    
    # Read coordinates
    coords = []
    for i in range(start_line, start_line + sum(atom_counts)):
        coords.append(list(map(float, lines[i].split()[:3])))
    coords = np.array(coords)
    
    # Check for velocities (optional in POSCAR)
    velocities = None
    if start_line + sum(atom_counts) < len(lines) and lines[start_line + sum(atom_counts)].strip():
        velocities = []
        for i in range(start_line + sum(atom_counts), start_line + 2 * sum(atom_counts)):
            velocities.append(list(map(float, lines[i].split()[:3])))
        velocities = np.array(velocities)
    
    return comment, lattice, elements, atom_counts, coords, coord_type, velocities

def fractional_to_cartesian(coords, lattice):
    """Convert fractional coordinates to Cartesian using lattice vectors."""
    return np.dot(coords, lattice)

def write_lammps_data(filename, comment, lattice, elements, atom_counts, coords, velocities):
    """Write LAMMPS data file in 'atom' format."""
    # Compute simulation box parameters
    a, b, c = lattice
    xlo, ylo, zlo = 0.0, 0.0, 0.0
    xhi = np.linalg.norm(a)
    yhi = np.linalg.norm(b)
    zhi = np.linalg.norm(c)
    
    # For non-orthorhombic cells, compute xy, xz, yz tilts
    xy = np.dot(b, a) / np.linalg.norm(a)
    xz = np.dot(c, a) / np.linalg.norm(a)
    yz = np.dot(c, b) / np.linalg.norm(b)
    
    # Assign atom types
    atom_types = []
    type_map = {elem: i+1 for i, elem in enumerate(set(elements))}
    for elem, count in zip(elements, atom_counts):
        atom_types.extend([type_map[elem]] * count)
    
    # Write LAMMPS data file
    with open(filename, 'w') as f:
        f.write(f"{comment} (Converted from POSCAR by poscar_to_lammps.py)\n\n")
        f.write(f"{sum(atom_counts)} atoms\n")
        f.write(f"{len(set(elements))} atom types\n\n")
        
        # Box dimensions
        f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
        f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
        f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n")
        if any([abs(xy) > 1e-6, abs(xz) > 1e-6, abs(yz) > 1e-6]):
            f.write(f"{xy:.6f} {xz:.6f} {yz:.6f} xy xz yz\n")
        f.write("\n")
        
        # Masses section
        f.write("Masses\n\n")
        for elem, type_id in sorted(type_map.items(), key=lambda x: x[1]):
            if elem not in ATOMIC_MASSES:
                raise ValueError(f"Atomic mass for element {elem} not found in ATOMIC_MASSES dictionary.")
            f.write(f"{type_id} {ATOMIC_MASSES[elem]:.6f} # {elem}\n")
        f.write("\n")
        
        # Atoms section
        f.write("Atoms # atom\n\n")
        for i, (type_id, coord) in enumerate(zip(atom_types, coords), 1):
            f.write(f"{i} {type_id} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")
        
        # Velocities section (if present)
        if velocities is not None:
            f.write("\nVelocities\n\n")
            for i, vel in enumerate(velocities, 1):
                f.write(f"{i} {vel[0]:.6f} {vel[1]:.6f} {vel[2]:.6f}\n")

def poscar_to_lammps(poscar_file, lammps_file):
    """Main function to convert POSCAR to LAMMPS data file."""
    comment, lattice, elements, atom_counts, coords, coord_type, velocities = read_poscar(poscar_file)
    
    # Convert coordinates if necessary
    if coord_type.startswith('d'):  # Direct (fractional) coordinates
        coords = fractional_to_cartesian(coords, lattice)
        if velocities is not None:
            velocities = fractional_to_cartesian(velocities, lattice)
    
    write_lammps_data(lammps_file, comment, lattice, elements, atom_counts, coords, velocities)
    print(f"LAMMPS data file written to {lammps_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python poscar_to_lammps.py <POSCAR_file> <LAMMPS_output_file>")
        sys.exit(1)
    poscar_to_lammps(sys.argv[1], sys.argv[2])