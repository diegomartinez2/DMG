# consolidate_forces.py
import numpy as np

# --- Parámetros de entrada ---
lammps_data_prefix = "displaced_config_"
lammps_outputs_dir = "lammps_outputs"
forces_prefix = "forces_"
dfset_output_file = "si_displacements.xyz"

# --- Funciones Auxiliares ---
def read_lammps_coords(filename):
    """Lee coordenadas de un archivo de datos de LAMMPS (formato 'Atoms' section)."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    atoms_section_start = -1
    for i, line in enumerate(lines):
        if line.strip() == "Atoms":
            atoms_section_start = i + 2 # +2 para saltar la línea "Atoms" y la línea en blanco
            break

    if atoms_section_start == -1:
        raise ValueError(f"Sección 'Atoms' no encontrada en {filename}")

    coords = []
    atom_types = []
    for line in lines[atoms_section_start:]:
        if not line.strip(): # Detener al final de la sección Atoms
            break
        parts = line.split()
        atom_types.append(parts[1]) # Guarda el tipo LAMMPS (1, 2, ...)
        coords.append([float(parts[2]), float(parts[3]), float(parts[4])]) # x, y, z

    return np.array(coords), atom_types

def read_lammps_forces(filename):
    """Lee fuerzas de un archivo de volcado de LAMMPS (dump custom id type x y z fx fy fz)."""
    # Saltamos las dos primeras líneas de encabezado de LAMMPS dump
    return np.loadtxt(filename, skiprows=2)[:, -3:] # Las últimas 3 columnas son fx, fy, fz

def map_lammps_type_to_element(lammps_type, unique_atoms_base):
    """Mapea el tipo de átomo de LAMMPS (1, 2, ...) a su símbolo de elemento (Si, O, etc.)."""
    # Esto depende de cómo se mapearon los tipos en write_lammps_data
    # Asume que el tipo 1 es el primer elemento único, tipo 2 el segundo, etc.
    unique_types = sorted(list(set(unique_atoms_base)))
    return unique_types[int(lammps_type) - 1]

# --- Lógica principal ---
print("Consolidando datos de fuerzas y desplazamientos...")

all_displacements_forces = []

# Asume que los números de configuración comienzan desde 1
i = 1
while True:
    config_num_str = f"{i:04d}"
    lammps_data_file = f"{lammps_data_prefix}{config_num_str}.lammps"
    forces_file = f"{lammps_outputs_dir}/{forces_prefix}{config_num_str}.dat"

    if not np.os.path.exists(lammps_data_file):
        if i == 1:
            print(f"Error: No se encontró el primer archivo de datos LAMMPS: {lammps_data_file}")
            print("Asegúrate de haber ejecutado 'generate_lammps_inputs.py' y 'run_all_lammps.sh'.")
            exit(1)
        break # Ya no hay más archivos
    if not np.os.path.exists(forces_file):
        print(f"Advertencia: No se encontraron fuerzas para {lammps_data_file}. Omitiendo.")
        i += 1
        continue

    coords, lammps_atom_types = read_lammps_coords(lammps_data_file)
    forces = read_lammps_forces(forces_file)

    # Necesitamos el mapeo de tipos LAMMPS (1,2..) a símbolos de elemento (Si, O..)
    # Para el ejemplo de Si, todos son Si, pero si hubiera varios tipos, lo necesitas del si.xyz original.
    # Por simplicidad, asumimos que todos los átomos son del mismo tipo que en si.xyz
    # Si tienes varios tipos, lee el si.xyz original y crea un mapeo de id_atom_lammps -> simbolo_elemento
    # Para este ejemplo simple de Si, podemos asumir el símbolo 'Si'
    original_atoms_base, _ = read_xyz("si.xyz") # Leer el original para obtener los símbolos
    unique_elements_base = sorted(list(set(original_atoms_base))) # ['Si'] en este caso

    all_displacements_forces.append({
        'coords': coords,
        'forces': forces,
        'comment': f"Configuration {config_num_str} from LAMMPS",
        'elements': [map_lammps_type_to_element(t, unique_elements_base) for t in lammps_atom_types]
    })
    i += 1

print(f"Consolidados {len(all_displacements_forces)} configuraciones.")

# --- Escribir el DFSET_FILE final ---
with open(dfset_output_file, 'w') as f:
    for config_data in all_displacements_forces:
        num_atoms = len(config_data['coords'])
        f.write(f"{num_atoms}\n")
        f.write(f"{config_data['comment']}\n")
        for j in range(num_atoms):
            elem = config_data['elements'][j]
            x, y, z = config_data['coords'][j]
            fx, fy, fz = config_data['forces'][j]
            f.write(f"{elem} {x:.6f} {y:.6f} {z:.6f} {fx:.6f} {fy:.6f} {fz:.6f}\n")

print(f"Archivo DFSET ({dfset_output_file}) generado exitosamente.")
