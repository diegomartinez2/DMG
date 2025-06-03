#!/usr/local/bin/python
import numpy as np

def read_poscar(poscar_file):
    with open(poscar_file, 'r') as f:
        lines = f.readlines()

    # Lese die Skalierungsfaktor
    scale = float(lines[1].strip())

    # Lese die Gittervektoren
    lattice_vectors = []
    for i in range(2, 5):
        lattice_vectors.append([float(x) for x in lines[i].split()])
    lattice_vectors = np.array(lattice_vectors) * scale

    # Lese die Atomtypen und deren Anzahl
    atom_types = lines[5].split()
    atom_counts = [int(x) for x in lines[6].split()]

    # Lese die Koordinaten
    coord_type = lines[7].strip().lower()
    coordinates = []
    for i in range(8, 8 + sum(atom_counts)):
        coordinates.append([float(x) for x in lines[i].split()[:3]])
    coordinates = np.array(coordinates)

    if coord_type == 'direct':
        # Umwandlung von direkten zu kartesischen Koordinaten
        coordinates = np.dot(coordinates, lattice_vectors)

    return lattice_vectors, atom_types, atom_counts, coordinates

def write_lammps_data(lattice_vectors, atom_types, atom_counts, coordinates, output_file):
    with open(output_file, 'w') as f:
        # Schreibe die Header-Informationen
        f.write("LAMMPS data file generated from POSCAR\n\n")
        f.write(f"{sum(atom_counts)} atoms\n")
        f.write(f"{len(atom_types)} atom types\n\n")

        # Schreibe die Box-Dimensionen
        f.write("0.0 {:.6f} xlo xhi\n".format(lattice_vectors[0, 0]))
        f.write("0.0 {:.6f} ylo yhi\n".format(lattice_vectors[1, 1]))
        f.write("0.0 {:.6f} zlo zhi\n".format(lattice_vectors[2, 2]))
        f.write("{:.6f} {:.6f} {:.6f} xy xz yz\n\n".format(lattice_vectors[1, 0], lattice_vectors[2, 0], lattice_vectors[2, 1]))

        # Schreibe die Atommassen (hier als Platzhalter)
        f.write("Masses\n\n")
        for i, atom_type in enumerate(atom_types):
            f.write(f"{i+1} 1.0  # {atom_type}\n")
        f.write("\n")

        # Schreibe die Atome
        f.write("Atoms\n\n")
        atom_id = 1
        for i, (atom_type, count) in enumerate(zip(atom_types, atom_counts)):
            for j in range(count):
                f.write(f"{atom_id} {i+1} {coordinates[atom_id-1, 0]:.6f} {coordinates[atom_id-1, 1]:.6f} {coordinates[atom_id-1, 2]:.6f}\n")
                atom_id += 1

# Beispielaufruf
poscar_file = 'POSCAR'
output
