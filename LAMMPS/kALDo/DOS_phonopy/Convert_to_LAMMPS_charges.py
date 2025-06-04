import os
import numpy as np

# Dictionary of charges for each atom type (modify as needed)
#CHARGE_DICT = {
#    "S": -2.0,   # Example: Sulfur charge
#    "C": 0.0,    # Example: Carbon charge
#    "Cu": 1.0    # Example: Copper charge
#}
CHARGE_DICT = {
    "S": -0.932707,   # Example: Sulfur charge
    "C": -0.0672927,    # Example: Carbon charge
    "Cu": 2.0    # Example: Copper charge
}

def read_lammps_data(filename):
    """Read a LAMMPS data file and extract box, atom types, masses, and coordinates."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Initialize variables
    natoms = 0
    ntypes = 0
    xlo, xhi, ylo, yhi, zlo, zhi = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    xy, xz, yz = 0.0, 0.0, 0.0
    masses = {}
    atom_types = []
    coords = []
    atom_type_ids = []

    # Parse header
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('#'):
            i += 1
            continue
        if 'atoms' in line:
            natoms = int(line.split()[0])
        if 'atom types' in line:
            ntypes = int(line.split()[0])
        if 'xlo xhi' in line:
            xlo, xhi = map(float, line.split()[:2])
        if 'ylo yhi' in line:
            ylo, yhi = map(float, line.split()[:2])
        if 'zlo zhi' in line:
            zlo, zhi = map(float, line.split()[:2])
        if 'xy xz yz' in line:
            xy, xz, yz = map(float, line.split()[:3])
        if line == 'Masses':
            i += 2  # Skip blank line
            for j in range(ntypes):
                type_id, mass, *comment = lines[i + j].split()
                type_id = int(type_id)
                mass = float(mass)
                # Extract element name from comment (e.g., '# S')
                element = comment[1] if len(comment) > 1 else f"Type{type_id}"
                masses[type_id] = (mass, element)
                atom_types.append(element)
            i += ntypes
        if line == 'Atoms # atomic':
            i += 2  # Skip blank line
            for j in range(natoms):
                atom_id, type_id, x, y, z = map(float, lines[i + j].split()[:5])
                coords.append([x, y, z])
                atom_type_ids.append(int(type_id))
            i += natoms
        i += 1

    # Convert coordinates to numpy array
    coords = np.array(coords)

    # Map atom type IDs to elements
    atom_list = [masses[type_id][1] for type_id in atom_type_ids]

    return (xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz), atom_types, masses, coords, atom_type_ids, atom_list

def write_lammps_data(filename, box, atom_types, masses, coords, atom_type_ids, atom_list):
    """Write LAMMPS data file with charges and zero velocities."""
    xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz = box
    natoms = len(coords)
    ntypes = len(atom_types)

    # Create atom type mapping (1-based indexing for LAMMPS)
    atom_type_map = {atype: i + 1 for i, atype in enumerate(atom_types)}

    with open(filename, 'w') as f:
        # Header
        f.write(f"LAMMPS data file generated from {filename}\n\n")
        f.write(f"{natoms} atoms\n")
        f.write(f"{ntypes} atom types\n\n")

        # Box dimensions
        f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
        f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
        f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n")
        if any([xy, xz, yz]):
            f.write(f"{xy:.6f} {xz:.6f} {yz:.6f} xy xz yz\n")
        f.write("\n")

        # Masses section
        f.write("Masses\n\n")
        for i, atype in enumerate(atom_types):
            type_id = i + 1
            mass = next(m for tid, (m, e) in masses.items() if e == atype)
            f.write(f"{type_id} {mass:.6f} # {atype}\n")
        f.write("\n")

        # Atoms section (atom_style charge: atom-ID, atom-type, charge, x, y, z)
        f.write("Atoms # charge\n\n")
        for i, (coord, atype, type_id) in enumerate(zip(coords, atom_list, atom_type_ids)):
            atom_type = atom_type_map[atype]
            charge = CHARGE_DICT.get(atype, 0.0)  # Default to 0.0 if not in dict
            f.write(f"{i + 1} {atom_type} {charge:.6f} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")
        f.write("\n")

        # Velocities section (all set to 0)
        f.write("Velocities\n\n")
        for i in range(natoms):
            f.write(f"{i + 1} 0.0 0.0 0.0\n")

def main():
    # Process files POSCAR-001.lmp to POSCAR-045.lmp
    for i in range(1, 46):
        input_file = f"POSCAR-{i:03d}.lmp"
        output_file = f"POSCAR-{i:03d}.lammps"

        if not os.path.exists(input_file):
            print(f"Warning: {input_file} not found, skipping.")
            continue

        try:
            box, atom_types, masses, coords, atom_type_ids, atom_list = read_lammps_data(input_file)
            write_lammps_data(output_file, box, atom_types, masses, coords, atom_type_ids, atom_list)
            print(f"Converted {input_file} to {output_file}")
        except Exception as e:
            print(f"Error processing {input_file}: {str(e)}")

if __name__ == "__main__":
    main()
