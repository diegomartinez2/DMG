#!/usr/bin/env python
import numpy as np
import os

# Constante de conversión de Angstroms a Bohr
# 1 Angstrom = 1.8897259886 Bohr
ANGSTROM_TO_BOHR = 1.8897259886

# --- Módulo 1: Parsing de información atómica y de red del SPOSCAR ---
def parse_sposcar_info(sposcar_filepath: str):
    """
    Parsea un archivo SPOSCAR para extraer información clave:
    número total de átomos, número de tipos de átomos, nombres de tipos de átomos únicos,
    la matriz de vectores de la red de la supercelda en Angstroms,
    la lista de todas las especies (como en línea 6), sus conteos (como en línea 7),
    las coordenadas atómicas y el tipo de coordenadas.

    Args:
        sposcar_filepath (str): Ruta al archivo SPOSCAR de la supercelda.

    Returns:
        tuple: (total_atoms, num_atom_kinds, unique_atom_species, supercell_lattice_matrix_angstrom,
                all_species_raw, all_counts_raw, atomic_coords_raw, coords_type_line),
               o (None, None, None, None, None, None, None, None) si hay un error.
    """
    try:
        with open(sposcar_filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo SPOSCAR en {sposcar_filepath}")
        return None, None, None, None, None, None, None, None
    except Exception as e:
        print(f"Error al leer el archivo SPOSCAR: {e}")
        return None, None, None, None, None, None, None, None

    try:
        # Línea 2: Factor de escala global (generalmente 1.0)
        global_scale_factor = float(lines[1].strip())

        # Líneas 3-5: Vectores de la red de la supercelda
        supercell_lattice_vectors_raw = []
        for i in range(2, 5):
            vec = list(map(float, lines[i].strip().split()))
            supercell_lattice_vectors_raw.append(vec)
        supercell_lattice_matrix_angstrom = np.array(supercell_lattice_vectors_raw) * global_scale_factor

        # Línea 6: Nombres de los tipos de átomos (pueden estar repetidos, ej. 'Si Si C')
        all_species_raw = lines[5].strip().split()

        # Eliminar duplicados manteniendo el orden para obtener los tipos únicos
        unique_atom_species = []
        [unique_atom_species.append(item) for item in all_species_raw if item not in unique_atom_species]

        # Línea 7: Conteo de átomos por tipo (correspondientes a all_species_raw, ej. '16 8 4')
        all_counts_raw = list(map(int, lines[6].strip().split()))

        # Calcular el número total de átomos en la supercelda
        total_atoms = sum(all_counts_raw)

        # El número de tipos de átomos es el número de especies únicas
        num_atom_kinds = len(unique_atom_species)

        # Línea 8: Tipo de coordenadas (Cartesian o Direct)
        coords_type_line = lines[7].strip()

        # Coordenadas atómicas (a partir de la línea 9 hasta el final de los átomos)
        atomic_coords_raw = []
        # La primera línea de coordenadas está en el índice 8 (0-indexed)
        # O sea, línea 9 del archivo.
        start_coord_line_idx = 8
        for i in range(start_coord_line_idx, start_coord_line_idx + total_atoms):
            coord = list(map(float, lines[i].strip().split()))
            atomic_coords_raw.append(coord)

        return total_atoms, num_atom_kinds, unique_atom_species, supercell_lattice_matrix_angstrom, \
               all_species_raw, all_counts_raw, atomic_coords_raw, coords_type_line

    except (IndexError, ValueError) as e:
        print(f"Error al parsear el SPOSCAR. Formato inesperado en las líneas: {e}")
        print("Asegúrate de que el SPOSCAR tiene el formato VASP tradicional (9 líneas de encabezado + coordenadas).")
        return None, None, None, None, None, None, None, None

# --- Módulo 2: Generación de la sección &general ---
def generate_general_section_content(total_atoms: int, num_atom_kinds: int, unique_atom_species: list, prefix: str = "displacements"):
    """
    Genera el contenido de la sección '&general'.

    Args:
        total_atoms (int): Número total de átomos.
        num_atom_kinds (int): Número de tipos de átomos.
        unique_atom_species (list): Lista de los nombres de los tipos de átomos únicos.
        prefix (str): Prefijo para los archivos de salida.

    Returns:
        str: Contenido formateado de la sección '&general'.
    """
    unique_atom_species_str = " ".join(unique_atom_species)
    content = ""
    content += "&general\n"
    content += f" PREFIX = {prefix} # ! Prefijo para los archivos de salida\n"
    content += f" MODE = suggest # ! Modo para generar patrones de desplazamiento\n"
    content += f" NAT = {total_atoms} # ! Número total de átomos en la supercelda\n"
    content += f" NKD = {num_atom_kinds} # ! Número de tipos de átomos\n"
    content += f" KD = {unique_atom_species_str} # ! Nombres de los tipos de átomos\n"
    content += "/\n"
    return content

# --- Módulo 3: Generación de la sección &interaction ---
def generate_interaction_section_content(norder: int = 2):
    """
    Genera el contenido de la sección '&interaction'.

    Args:
        norder (int): Orden de las interacciones anarmónicas a incluir (1: armónico, 2: cúbico).

    Returns:
        str: Contenido formateado de la sección '&interaction'.
    """
    content = ""
    content += "&interaction\n"
    content += f" NORDER = {norder} # 1: armónico, 2: cúbico\n"
    content += "/\n"
    return content

# --- Módulo 4: Generación de la sección &cell ---
def generate_cell_section_content(supercell_lattice_matrix_angstrom: np.ndarray):
    """
    Genera el contenido de la sección '&cell'.

    Args:
        supercell_lattice_matrix_angstrom (np.ndarray): Matriz de 3x3 de los vectores de la red
                                                        de la supercelda en Angstroms.

    Returns:
        str: Contenido formateado de la sección '&cell'.
    """
    # Calcular el factor de escala de ALAMODE en Bohr
    alamode_scale_factor_angstrom = np.linalg.norm(supercell_lattice_matrix_angstrom[0])
    alamode_scale_factor_bohr = alamode_scale_factor_angstrom * ANGSTROM_TO_BOHR

    # Calcular los vectores de la red normalizados para ALAMODE
    normalized_lattice_vectors = supercell_lattice_matrix_angstrom / alamode_scale_factor_bohr

    content = ""
    content += "&cell\n"
    content += "# !--------------------------------------------------------------------------------------------------\n"
    content += f" {alamode_scale_factor_bohr:.16f} # ! Factor de escala de la red en Bohr (Constante de red de la supercelda)\n"
    for vec in normalized_lattice_vectors:
        content += f" {vec[0]:.16f} {vec[1]:.16f} {vec[2]:.16f}\n"
    content += "/\n"
    return content

# --- Módulo 5: Generación de la sección &cutoff ---
def generate_cutoff_section_content(cutoff_2body: float = 5.0, cutoff_3body: float = 3.0):
    """
    Genera el contenido de la sección '&cutoff' usando el formato simplificado de asteriscos.

    Args:
        cutoff_2body (float): Radio de corte para interacciones de dos cuerpos.
        cutoff_3body (float): Radio de corte para interacciones de tres cuerpos.

    Returns:
        str: Contenido formateado de la sección '&cutoff'.
    """
    content = ""
    content += "&cutoff\n"
    content += "# !--------------------------------------------------------------------------------------------------\n"
    content += f" *-* {cutoff_2body:.1f} {cutoff_3body:.1f}\n" # Formato simplificado
    content += "/\n"
    return content

# --- Módulo 6: Generación de la sección &position ---
def generate_position_section_content(atomic_coords_raw: list, unique_atom_species: list,
                                      all_species_raw: list, all_counts_raw: list,
                                      coords_type_line: str):
    """
    Genera el contenido de la sección '&position'.

    Args:
        atomic_coords_raw (list): Lista de listas con las coordenadas [x, y, z] de cada átomo.
        unique_atom_species (list): Lista de los nombres de los tipos de átomos únicos (ej. ['Si', 'C']).
        all_species_raw (list): Lista de los nombres de los tipos de átomos tal como aparecen en la línea 6 del SPOSCAR.
        all_counts_raw (list): Lista de los conteos de átomos tal como aparecen en la línea 7 del SPOSCAR.
        coords_type_line (str): La línea del SPOSCAR que indica el tipo de coordenadas ('Cartesian' o 'Direct').

    Returns:
        str: Contenido formateado de la sección '&position'.
    """
    content = ""
    content += "&position\n"

    # Advertencia si las coordenadas no son cartesianas
    if coords_type_line.strip().lower() != "cartesian":
        content += f"# ! ADVERTENCIA: Las coordenadas del SPOSCAR son '{coords_type_line}'. ALAMODE espera coordenadas cartesianas.\n"
        content += "# ! Asegúrate de que estas coordenadas son interpretables por ALAMODE como cartesianas o conviértelas previamente.\n"

    # Crear un mapeo de nombre de especie a ID (1-basado)
    atom_id_map = {species: i + 1 for i, species in enumerate(unique_atom_species)}

    # Construir la secuencia de IDs de tipo de átomo para cada átomo en el orden del SPOSCAR
    atom_type_sequence = []
    current_atom_idx = 0
    for i, species_name_block in enumerate(all_species_raw):
        count_in_block = all_counts_raw[i]
        atom_type_id = atom_id_map[species_name_block]
        atom_type_sequence.extend([atom_type_id] * count_in_block)

    # Asegurarse de que el número de coordenadas coincide con la secuencia de tipos
    if len(atomic_coords_raw) != len(atom_type_sequence):
        print("Error: El número de coordenadas atómicas no coincide con el número total de átomos inferido de los conteos.")
        return "" # Retornar vacío o manejar el error de otra forma

    # Escribir cada línea de posición
    for i, coord in enumerate(atomic_coords_raw):
        type_id = atom_type_sequence[i]
        content += f" {type_id} {coord[0]:.16f} {coord[1]:.16f} {coord[2]:.16f}\n"
    content += "/\n"
    return content

# --- Función principal para generar el archivo alm_suggest.in ---
def generate_alm_suggest_in(sposcar_filepath: str, output_filename: str = "alm_suggest.in", prefix: str = "displacements"):
    """
    Genera el archivo 'alm_suggest.in' directamente a partir de un archivo SPOSCAR.

    Args:
        sposcar_filepath (str): Ruta al archivo SPOSCAR de la supercelda.
        output_filename (str): Nombre del archivo de salida a generar.
        prefix (str): Prefijo para los archivos de salida de ALAMODE.
    """
    # Parsear toda la información relevante del SPOSCAR
    total_atoms, num_atom_kinds, unique_atom_species, supercell_lattice_matrix_angstrom, \
    all_species_raw, all_counts_raw, atomic_coords_raw, coords_type_line = \
        parse_sposcar_info(sposcar_filepath)

    if total_atoms is None:
        return # Salir si hubo un error al parsear

    # Generar contenido de las secciones
    general_section = generate_general_section_content(total_atoms, num_atom_kinds, unique_atom_species, prefix)
    interaction_section = generate_interaction_section_content()
    cell_section = generate_cell_section_content(supercell_lattice_matrix_angstrom)
    cutoff_section = generate_cutoff_section_content()
    position_section = generate_position_section_content(atomic_coords_raw, unique_atom_species,
                                                         all_species_raw, all_counts_raw,
                                                         coords_type_line)

    # Escribir el archivo
    try:
        with open(output_filename, 'w') as out_f:
            out_f.write(general_section)
            out_f.write("\n")
            out_f.write(interaction_section)
            out_f.write("\n")
            out_f.write(cell_section)
            out_f.write("\n")
            out_f.write(cutoff_section)
            out_f.write("\n")
            out_f.write(position_section)
            out_f.write("\n")

            # Marcador para futuras secciones
            out_f.write("! --- Aquí irán las secciones &species y &atom --- \n")
            out_f.write("! Aún faltan estas secciones para completar el archivo.\n")


        print(f"\n¡Éxito! Archivo '{output_filename}' generado con las secciones '&general', '&interaction', '&cell', '&cutoff' y '&position'.")
        print("Por favor, revisa el contenido de este archivo antes de usarlo con ALAMODE.")

    except Exception as e:
        print(f"Error al escribir el archivo '{output_filename}': {e}")

# --- Bloque principal de ejecución ---
if __name__ == "__main__":
    print("--- Generador Modular de 'alm_suggest.in' (Secciones &general, &interaction, &cell, &cutoff, &position) ---")

    sposcar_input_file = input("Introduce la ruta al archivo SPOSCAR de la supercelda: ")

    # Puedes cambiar el prefijo si lo deseas
    output_prefix = "my_alamode_run"

    generate_alm_suggest_in(sposcar_input_file, prefix=output_prefix)
