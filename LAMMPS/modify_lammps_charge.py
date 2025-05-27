import uuid

# Define charges for each atom type (hypothetical values)
charges = {
    1: -0.932707,  # S (Sulfur)
    2: 0.0672927,   # C (Carbon)
    3: 2    # Cu (Copper)
}

def read_lammps_file(input_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    return lines

def write_lammps_file(output_file, header, masses, atoms):
    with open(output_file, 'w') as f:
        # Write header, replacing 'atomic' with 'charge' in the Atoms section
        for line in header:
            if line.strip().startswith("Atoms"):
                f.write("Atoms # charge\n")
            else:
                f.write(line)
        f.write("\n")
        # Write masses
        for line in masses:
            f.write(line)
        f.write("\n")
        # Write atoms with charge
        f.write("Atoms # charge\n\n")
        for atom in atoms:
            atom_id, atom_type, x, y, z = atom
            charge = charges.get(int(atom_type), 0.0)  # Default to 0.0 if type not found
            # Convert x, y, z to floats
            #x, y, z = float(x), float(y), float(z)
            try:
                x, y, z = float(x), float(y), float(z)
            except ValueError as e:
                print(f"Error converting coordinates for atom {atom_id}: {e}")
                raise
            f.write(f"{atom_id:>10} {atom_type:>4} {charge:>10.6f} {x:>12.6f} {y:>12.6f} {z:>12.6f} 0 0 0\n")
        f.write("Velocities\n\n")
        for atom in atoms:
            f.write(f"{atom_id:>10} 0.0 0.0 0.0\n")

def parse_lammps_file(lines):
    header = []
    masses = []
    atoms = []
    section = None

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            if section is None:
                header.append(line + "\n")
            continue

        if line.startswith("Masses"):
            section = "masses"
            masses.append(line + "\n")
            continue
        elif line.startswith("Atoms"):
            section = "atoms"
            continue

        if section == "masses" and line:
            masses.append(line + "\n")
        elif section == "atoms" and line:
            parts = line.split()
            if len(parts) >= 5:
                atom_id, atom_type, x, y, z = parts[:5]
                atoms.append((atom_id, atom_type, x, y, z))
        elif section is None:
            header.append(line + "\n")

    return header, masses, atoms

def main(input_file, output_file):
    # Read and parse the input file
    lines = read_lammps_file(input_file)
    header, masses, atoms = parse_lammps_file(lines)

    # Write the modified file
    write_lammps_file(output_file, header, masses, atoms)
    print(f"Modified LAMMPS file written to {output_file}")

if __name__ == "__main__":
    input_file = "input.lammps"
    output_file = "output.lammps"
    main(input_file, output_file)
