#!/usr/bin/env python
# generate_lammps_inputs.py
import numpy as np

# --- Parámetros de entrada ---
xyz_file = "si.xyz"
harmonic_pattern_file = "Si_patterns.HARMONIC_pattern"
cubic_pattern_file = "Si_patterns.ANHARM2_pattern" # O ANHARM3_pattern si NORDER=3
displacement_magnitude = 0.01 # Magnitud del desplazamiento en Angstroms (ej. 0.01 Å)
lammps_data_prefix = "displaced_config_"
potential_name = "Si" # Nombre de tu potencial para el tipo de atomo en LAMMPS (ej. Si para Stillinger-Weber)

# --- Funciones Auxiliares ---
def read_xyz(filename):
    """Lee un archivo .xyz y devuelve átomos y coordenadas."""
    with open(filename, 'r') as f:
        num_atoms = int(f.readline())
        comment = f.readline().strip()
        atoms = []
        coords = []
        for _ in range(num_atoms):
            line = f.readline().split()
            atoms.append(line[0])
            coords.append([float(x) for x in line[1:4]])
    return atoms, np.array(coords)

def write_lammps_data(filename, atoms, coords, box_dims):
    """Escribe un archivo de datos de LAMMPS."""
    with open(filename, 'w') as f:
        f.write(f"LAMMPS data file - Displaced configuration\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{len(set(atoms))} atom types\n\n") # Asume que los tipos son 1, 2, ...

        f.write(f"0.0 {box_dims[0]:.6f} xlo xhi\n")
        f.write(f"0.0 {box_dims[1]:.6f} ylo yhi\n")
        f.write(f"0.0 {box_dims[2]:.6f} zlo zhi\n\n")

        f.write("Masses\n\n")
        unique_atoms = sorted(list(set(atoms)))
        for i, atom_type in enumerate(unique_atoms):
            # Necesitarías la masa real del átomo aquí, ej. Si = 28.0855
            # Para este ejemplo, solo usamos un marcador de posición
            f.write(f"{i+1} 28.0855 # {atom_type}\n")

        f.write("\nAtoms\n\n")
        atom_type_map = {atom_type: i+1 for i, atom_type in enumerate(unique_atoms)}
        for i, (atom_sym, coord) in enumerate(zip(atoms, coords)):
            f.write(f"{i+1} {atom_type_map[atom_sym]} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")

def read_pattern_file(filename):
    """Lee un archivo .pattern de ALAMODE."""
    patterns = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            num_displacements = int(lines[i].strip())
            i += 1
            displacements = []
            for _ in range(num_displacements):
                parts = lines[i].split()
                # Atom_id (1-based), dx, dy, dz
                displacements.append((int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])))
                i += 1
            patterns.append(displacements)
    return patterns

# --- Lógica principal ---
print(f"Leyendo estructura base de {xyz_file}...")
atoms_base, coords_base = read_xyz(xyz_file)

# Para un archivo .xyz, asumimos que es una celda unitaria.
# Necesitas construir una supercelda para LAMMPS.
# Para el ejemplo de Si, se suele usar una supercelda 2x2x2 o 3x3x3.
# Este script NO construye la supercelda, solo desplaza los átomos base.
# Si tu si.xyz ya es una supercelda, está bien.
# Si no, necesitarías expandir la celda unitaria a una supercelda aquí.
# Para simplificar, asumiremos que si.xyz ya representa la supercelda para los desplazamientos.

# Estimación de las dimensiones de la caja (asume una caja cúbica para Si si.xyz simple)
# Si tu si.xyz es una celda unitaria, tendrás que ajustarlo
# Si es una supercelda, puedes estimar con el rango de coordenadas
box_x = np.max(coords_base[:,0]) - np.min(coords_base[:,0])
box_y = np.max(coords_base[:,1]) - np.min(coords_base[:,1])
box_z = np.max(coords_base[:,2]) - np.min(coords_base[:,2])
# En el caso de Si (cubico), puedes usar el parámetro de red real para la supercelda
# Ejemplo: para Si, a = 5.43 Angstroms. Para una supercelda 2x2x2, seria 2*5.43.
# Si no tienes esa información aquí, usa las estimaciones de max/min como base.
# ESTO ES CRÍTICO: Las dimensiones de la caja deben ser precisas para LAMMPS.
# Si conoces el parámetro de red y la multiplicidad de la supercelda:
# por ejemplo, para Si celda de diamante, a=5.43. Una supercelda 2x2x2 tendría 16 átomos y dim=10.86
# box_x = box_y = box_z = 10.86 # Ejemplo para Si 2x2x2 supercelda

box_dims = [box_x, box_y, box_z] # Ajusta esto a las dimensiones REALES de tu supercelda


pattern_files = [harmonic_pattern_file, cubic_pattern_file]
all_patterns = []
for p_file in pattern_files:
    if np.os.path.exists(p_file):
        all_patterns.extend(read_pattern_file(p_file))
        print(f"Cargado {len(read_pattern_file(p_file))} patrones de {p_file}")
    else:
        print(f"Advertencia: Archivo de patrón no encontrado: {p_file}. Omitiendo.")

print(f"Generando {len(all_patterns)} configuraciones desplazadas para LAMMPS...")
for i, pattern in enumerate(all_patterns):
    coords_displaced = np.copy(coords_base)
    for atom_id, dx, dy, dz in pattern:
        # ALM usa indices 1-based, Python usa 0-based
        coords_displaced[atom_id - 1] += np.array([dx, dy, dz]) * displacement_magnitude

    output_filename = f"{lammps_data_prefix}{i+1:04d}.lammps"
    write_lammps_data(output_filename, atoms_base, coords_displaced, box_dims)
    print(f"Creado {output_filename}")

print("Generación de archivos de entrada de LAMMPS completada.")
