#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Construir_aml_anphon_file.py
#
#  Copyright 2025 Diego Martinez Gutierrez
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
'''
Este programa crea los input files para ALAMODE.
'''
# ---------------------------
# Importación de los módulos
# ---------------------------

import numpy as np
import os

# ----------
# Parametros
# ----------
# Constante de conversión de Angstroms a Bohr
# 1 Angstrom = 1.8897259886 Bohr
ANGSTROM_TO_BOHR = 1.8897259886

# Diccionario de masas atómicas (ejemplos, expandir según sea necesario)
# Fuente: Pesos atómicos estándar de la IUPAC (aproximados para la mayoría de los casos)
ATOMIC_MASSES = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305,
    "Al": 26.982, "Si": 28.0855, "P": 30.974, "S": 32.06, "Cl": 35.45, "Ar": 39.948,
    "K": 39.0983, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996,
    "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.63, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.96,
    "Tc": 98.0, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41,
    "In": 114.82, "Sn": 118.71, "Sb": 121.76, "I": 126.90, "Te": 127.60, "Xe": 131.29,
    "Cs": 132.91, "Ba": 137.33, "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24,
    "Pm": 145.0, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.50,
    "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97, "Hf": 178.49,
    "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08,
    "Au": 196.97, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Po": 209.0,
    "At": 210.0, "Rn": 222.0, "Fr": 223.0, "Ra": 226.0, "Ac": 227.0, "Th": 232.04,
    "Pa": 231.04, "U": 238.03, "Np": 237.0, "Pu": 244.0, "Am": 243.0, "Cm": 247.0,
    "Bk": 247.0, "Cf": 251.0, "Es": 252.0, "Fm": 257.0, "Md": 258.0, "No": 259.0,
    "Lr": 262.0, "Rf": 267.0, "Db": 268.0, "Sg": 271.0, "Bh": 272.0, "Hs": 270.0,
    "Mt": 276.0, "Ds": 281.0, "Rg": 280.0, "Cn": 285.0, "Nh": 286.0, "Fl": 289.0,
    "Mc": 290.0, "Lv": 293.0, "Ts": 294.0, "Og": 294.0
}

# ----------
# Funciones
# ----------

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
def generate_general_section_content(total_atoms: int, num_atom_kinds: int, unique_atom_species: list, prefix: str = "run", program_type: str = "alm", program_mode: str = "suggest", fcsxml_file: str = None):
    """
    Genera el contenido de la sección '&general'.

    Args:
        total_atoms (int): Número total de átomos.
        num_atom_kinds (int): Número de tipos de átomos.
        unique_atom_species (list): Lista de los nombres de los tipos de átomos únicos.
        prefix (str): Prefijo para los archivos de salida.
        program_type (str): El programa objetivo ('alm' o 'anphon').
        program_mode (str): El modo de operación específico del programa (ej. 'suggest', 'optimize' para ALAMODE; 'phonons', 'RTA', 'SCPH' para ANPHON).
        fcsxml_file (str, optional): Nombre del archivo XML de constantes de fuerza (FCSXML). Necesario para ANPHON.

    Returns:
        str: Contenido formateado de la sección '&general'.
    """
    unique_atom_species_str = " ".join(unique_atom_species)
    content = ""
    content += "&general\n"
    content += f" PREFIX = {prefix} # ! Prefijo para los archivos de salida\n"
    content += f" MODE = {program_mode} # ! Modo de operación de {program_type.upper()}\n"

    if program_type.lower() == "anphon":
        if fcsxml_file:
            content += f" FCSXML = {fcsxml_file} # ! Archivo XML de constantes de fuerza\n"
        else:
            content += "# ! ADVERTENCIA: FCSXML es esencial para ANPHON. Asegúrate de que existe o añádelo manualmente.\n"

        # Añadir la línea MASS para ANPHON
        mass_values = []
        for species in unique_atom_species:
            mass = ATOMIC_MASSES.get(species)
            if mass is None:
                print(f"¡ADVERTENCIA! Masa para el elemento '{species}' no encontrada en el diccionario. Por favor, añádela manualmente.")
                mass_values.append("???") # Placeholder para que el usuario sepa que falta
            else:
                mass_values.append(f"{mass:.4f}")
        content += f" MASS = {' '.join(mass_values)} # ! Masas atómicas de los elementos definidos en KD\n"


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
    # Calcular el factor de escala de ALAMODE/ANPHON en Bohr
    alamode_scale_factor_angstrom = np.linalg.norm(supercell_lattice_matrix_angstrom[0])
    alamode_scale_factor_bohr = alamode_scale_factor_angstrom * ANGSTROM_TO_BOHR

    # Calcular los vectores de la red normalizados para ALAMODE/ANPHON
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

# --- Módulo 7: Generación de la sección &optimize ---
def generate_optimize_section_content(optimize_type: str = "harmonic", fc2xml_file: str = None):
    """
    Genera el contenido de la sección '&optimize'.

    Args:
        optimize_type (str): Tipo de optimización ('harmonic' para DFSET_harmonic, 'cubic' para DFSET_cubic).
                              Por defecto: 'harmonic'.
        fc2xml_file (str, optional): Nombre del archivo XML de IFCs de segundo orden si
                                     se usa 'DFSET_cubic'. Por defecto: None.

    Returns:
        str: Contenido formateado de la sección '&optimize'.
    """
    content = ""
    content += "&optimize\n"

    if optimize_type.lower() == "harmonic":
        content += f" DFSET = DFSET_harmonic\n"
    elif optimize_type.lower() == "cubic":
        content += f" DFSET = DFSET_cubic\n"
        if fc2xml_file:
            content += f" FC2XML = {fc2_xml_file} # Fix harmonic IFCs\n"
        else:
            content += f"# ! ADVERTENCIA: Para DFSET_cubic, se recomienda especificar FC2XML.\n"
    else:
        print(f"Advertencia: Tipo de optimización '{optimize_type}' no reconocido. Usando DFSET_harmonic por defecto.")
        content += f" DFSET = DFSET_harmonic\n"

    content += "/\n"
    return content

# --- Módulo 8: Generación de la sección &kpoint ---
def generate_kpoint_section_content(kpmode: int = 2, kpoint_params: list = None):
    """
    Genera el contenido de la sección '&kpoint' para ANPHON.

    Args:
        kpmode (int): Modo de puntos k (1: line mode, 2: uniform mesh mode). Por defecto: 2.
        kpoint_params (list): Parámetros específicos del modo.
                              Si kpmode = 1: Lista de strings de líneas de banda (ej. ["G 0.0 0.0 0.0 X 0.5 0.5 0.0 51"]).
                              Si kpmode = 2: Lista de 3 enteros [nx, ny, nz] para la malla uniforme.
                              Por defecto: None, que resultará en [20, 20, 20] para kpmode 2.

    Returns:
        str: Contenido formateado de la sección '&kpoint'.
    """
    content = ""
    content += "&kpoint\n"
    content += f" {kpmode} # KPMODE = {kpmode}: {'line mode' if kpmode == 1 else 'uniform mesh mode'}\n"

    if kpmode == 1: # Line mode
        if not kpoint_params:
            # Ejemplo por defecto si no se dan parámetros para modo línea (ejemplo del manual para Si)
            content += " G 0.0 0.0 0.0 X 0.5 0.5 0.0 51\n"
            content += " X 0.5 0.5 1.0 G 0.0 0.0 0.0 51\n"
            content += " G 0.0 0.0 0.0 L 0.5 0.5 0.5 51\n"
            content += "# ! ADVERTENCIA: Se usaron líneas de banda por defecto. Modificar según necesidad.\n"
        else:
            for line_def in kpoint_params:
                content += f" {line_def}\n"
    elif kpmode == 2: # Uniform mesh mode
        if not kpoint_params or len(kpoint_params) != 3:
            default_mesh = [20, 20, 20]
            content += f" {default_mesh[0]} {default_mesh[1]} {default_mesh[2]}\n"
            content += "# ! ADVERTENCIA: Se usó una malla uniforme por defecto (20 20 20). Modificar según necesidad.\n"
        else:
            content += f" {kpoint_params[0]} {kpoint_params[1]} {kpoint_params[2]}\n"
    else:
        print(f"Advertencia: KPMODE = {kpmode} no es un modo válido. Usando 2 (uniform mesh) con valores por defecto.")
        content += f" 2 # KPMODE = 2: uniform mesh mode\n"
        content += f" 20 20 20\n"
        content += "# ! ADVERTENCIA: Se usó una malla uniforme por defecto (20 20 20) debido a KPMODE inválido.\n"

    content += "/\n"
    return content


# --- Función principal para generar el archivo alm_*.in / anphon_*.in ---
def generate_input_file(sposcar_filepath: str, output_filename: str, prefix: str = "run", program_type: str = "alm", program_mode: str = "suggest", optimize_type: str = "harmonic", fc2xml_file: str = None, fcsxml_file: str = None, kpmode: int = 2, kpoint_params: list = None):
    """
    Genera el archivo de entrada para ALAMODE o ANPHON directamente a partir de un archivo SPOSCAR.

    Args:
        sposcar_filepath (str): Ruta al archivo SPOSCAR de la supercelda.
        output_filename (str): Nombre del archivo de salida a generar (ej. alm_suggest.in, anphon_phonons.in).
        prefix (str): Prefijo para los archivos de salida de ALAMODE/ANPHON.
        program_type (str): El programa objetivo ('alm' o 'anphon').
        program_mode (str): Modo de operación específico del programa.
                            Para ALAMODE: 'suggest' o 'optimize'.
                            Para ANPHON: 'phonons', 'RTA', 'SCPH'.
        optimize_type (str): Tipo de optimización para ALAMODE 'optimize' ('harmonic' o 'cubic').
        fc2xml_file (str, optional): Nombre del archivo XML de IFCs de segundo orden para ALAMODE 'optimize' con cubic.
        fcsxml_file (str, optional): Nombre del archivo XML de constantes de fuerza (FCSXML) para ANPHON.
        kpmode (int): Modo de puntos k para ANPHON (1: line mode, 2: uniform mesh mode).
        kpoint_params (list): Parámetros específicos de puntos k para ANPHON.
    """
    # Parsear toda la información relevante del SPOSCAR
    total_atoms, num_atom_kinds, unique_atom_species, supercell_lattice_matrix_angstrom, \
    all_species_raw, all_counts_raw, atomic_coords_raw, coords_type_line = \
        parse_sposcar_info(sposcar_filepath)

    if total_atoms is None:
        return # Salir si hubo un error al parsear

    # --- Generar contenido de las secciones comunes ---
    general_section = generate_general_section_content(total_atoms, num_atom_kinds, unique_atom_species, prefix, program_type, program_mode, fcsxml_file)
    cell_section = generate_cell_section_content(supercell_lattice_matrix_angstrom)
    position_section = generate_position_section_content(atomic_coords_raw, unique_atom_species,
                                                         all_species_raw, all_counts_raw,
                                                         coords_type_line)

    # --- Generar contenido de secciones específicas por programa ---
    interaction_section = ""
    cutoff_section = ""
    optimize_section = ""
    kpoint_section = ""

    if program_type.lower() == "alm":
        interaction_section = generate_interaction_section_content()
        cutoff_section = generate_cutoff_section_content()
        if program_mode.lower() == "optimize":
            optimize_section = generate_optimize_section_content(optimize_type, fc2xml_file)
    elif program_type.lower() == "anphon":
        kpoint_section = generate_kpoint_section_content(kpmode, kpoint_params)
        # ANPHON también puede necesitar &interaction y &cutoff si usamos anarmonicidades (DFSET_cubic/quartic)
        # Por ahora, solo añadimos kpoint y FCSXML en general.
        pass

    # Escribir el archivo
    try:
        with open(output_filename, 'w') as out_f:
            out_f.write(general_section)
            out_f.write("\n")
            out_f.write(cell_section)
            out_f.write("\n")

            # Secciones específicas de ALAMODE
            if interaction_section: out_f.write(interaction_section + "\n")
            if cutoff_section: out_f.write(cutoff_section + "\n")

            out_f.write(position_section)
            out_f.write("\n")

            if optimize_section: out_f.write(optimize_section + "\n")

            # Secciones específicas de ANPHON
            if kpoint_section: out_f.write(kpoint_section + "\n")

            # Marcador para futuras secciones
            out_f.write("! --- Aún faltan las secciones &species y &atom (¡próximamente!) --- \n")
            out_f.write("! Si es ANPHON, también faltan: &scph, &qha, &relax, &strain, &displace.\n")


        print(f"\n¡Éxito! Archivo '{output_filename}' generado.")
        print("Por favor, revisa el contenido de este archivo antes de usarlo con ALAMODE/ANPHON.")

    except Exception as e:
        print(f"Error al escribir el archivo '{output_filename}': {e}")

# --- Bloque principal de ejecución ---
if __name__ == "__main__":
    print("--- Generador Modular de archivos de entrada para ALAMODE/ANPHON ---")

    sposcar_input_file = input("Introduce la ruta al archivo SPOSCAR de la supercelda: ")

    # --- Elección del programa ---
    program_choice = input("¿Deseas generar el archivo para 'alm' o 'anphon'? (alm/anphon): ").strip().lower()
    if program_choice not in ["alm", "anphon"]:
        print("Programa no válido. Se usará 'alm' por defecto.")
        program_choice = "alm"

    # --- Variables para modos y archivos XML ---
    alm_mode = "suggest"
    optimize_calculation_type = "harmonic"
    fc2_xml_file = None

    anphon_mode = "phonons"
    anphon_fcs_xml_file = None
    anphon_kpmode = 2
    anphon_kpoint_params = None

    # --- Configuración detallada según el programa ---
    if program_choice == "alm":
        alm_mode_choice = input("¿Modo de ALAMODE? (suggest/optimize): ").strip().lower()
        if alm_mode_choice in ["suggest", "optimize"]:
            alm_mode = alm_mode_choice
        else:
            print("Modo ALAMODE no válido. Se usará 'suggest' por defecto.")

        if alm_mode == "optimize":
            opt_type_choice = input("¿Tipo de optimización: 'harmonic' o 'cubic'? (harmonic/cubic): ").strip().lower()
            if opt_type_choice == "cubic":
                optimize_calculation_type = "cubic"
                fc2_xml_file = input("Introduce el nombre del archivo FC2XML (ej. si222.xml) o deja en blanco si no lo tienes: ").strip()
                if not fc2_xml_file:
                    fc2_xml_file = None
            elif opt_type_choice != "harmonic":
                print("Tipo de optimización no válido. Se usará 'harmonic' por defecto.")
                optimize_calculation_type = "harmonic"

    elif program_choice == "anphon":
        anphon_mode_choice = input("¿Modo de ANPHON? (phonons/RTA/SCPH): ").strip().lower()
        if anphon_mode_choice in ["phonons", "rta", "scph"]:
            anphon_mode = anphon_mode_choice
        else:
            print("Modo ANPHON no válido. Se usará 'phonons' por defecto.")
            anphon_mode = "phonons"

        anphon_fcs_xml_file = input("Introduce el nombre del archivo FCSXML (ej. si222.xml) para ANPHON: ").strip()
        if not anphon_fcs_xml_file:
            print("¡ADVERTENCIA! FCSXML es crucial para ANPHON. Se dejará vacío por ahora.")
            anphon_fcs_xml_file = None

        kpmode_choice = input("¿Modo de puntos k para ANPHON? (1: línea de banda, 2: malla uniforme): ").strip()
        try:
            anphon_kpmode = int(kpmode_choice)
            if anphon_kpmode not in [1, 2]:
                raise ValueError
        except ValueError:
            print("KPMODE no válido. Se usará 2 (malla uniforme) por defecto.")
            anphon_kpmode = 2

        if anphon_kpmode == 1:
            print("\nIntroduce las líneas de banda (una por línea, ej. 'G 0.0 0.0 0.0 X 0.5 0.5 0.0 51').")
            print("Deja una línea en blanco y presiona Enter para finalizar la entrada.")
            lines = []
            while True:
                line = input(f"Línea de banda {len(lines)+1}: ").strip()
                if not line:
                    break
                lines.append(line)
            anphon_kpoint_params = lines if lines else None
            if not anphon_kpoint_params:
                print("No se introdujeron líneas de banda. Se usarán valores por defecto para KPMODE 1.")
        elif anphon_kpmode == 2:
            mesh_input = input("Introduce los tres enteros para la malla uniforme (ej. '20 20 20') o deja en blanco para 20 20 20: ").strip().split()
            if len(mesh_input) == 3:
                try:
                    anphon_kpoint_params = [int(n) for n in mesh_input]
                except ValueError:
                    print("Entrada de malla no válida. Se usará 20 20 20 por defecto.")
                    anphon_kpoint_params = None
            else:
                print("Entrada de malla no válida. Se usará 20 20 20 por defecto.")
                anphon_kpoint_params = None

    # --- Determinación del nombre de archivo de salida y prefijo ---
    output_prefix = "my_run" # Prefijo genérico, se puede personalizar

    if program_choice == "alm":
        output_file_name = f"alm_{alm_mode}.in"
        final_program_mode = alm_mode
    elif program_choice == "anphon":
        output_file_name = f"anphon_{anphon_mode}.in"
        final_program_mode = anphon_mode
    else:
        output_file_name = "unknown_program.in"
        final_program_mode = "default"
        print("Error interno: Tipo de programa desconocido.")


    # --- Llamada a la función principal ---
    generate_input_file(sposcar_input_file,
                            output_file_name,
                            prefix=output_prefix,
                            program_type=program_choice,
                            program_mode=final_program_mode,
                            optimize_type=optimize_calculation_type,
                            fc2xml_file=fc2_xml_file,
                            fcsxml_file=anphon_fcs_xml_file,
                            kpmode=anphon_kpmode,
                            kpoint_params=anphon_kpoint_params)
