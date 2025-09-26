#!/usr/bin/env python
import numpy as np
import os

# Constante de conversión de Angstroms a Bohr
# 1 Angstrom = 1.8897259886 Bohr
ANGSTROM_TO_BOHR = 1.8897259886

def get_alamode_cell_parameters(poscar_filepath: str, supercell_dim: tuple = (2, 2, 2)):
    """
    Calcula el factor de escala en Bohr y los vectores de red normalizados
    para el archivo de entrada de ALAMODE (alm_suggest.in) a partir de un POSCAR
    de la celda unitaria y una dimensión de supercelda.

    Args:
        poscar_filepath (str): Ruta al archivo POSCAR de la celda unitaria.
        supercell_dim (tuple): Dimensiones de la supercelda (e.g., (2, 2, 2) para 2x2x2).

    Returns:
        tuple: (alamode_scale_factor_bohr, normalized_lattice_vectors),
               o (None, None) si hay un error.
    """
    try:
        with open(poscar_filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo POSCAR en {poscar_filepath}")
        return None, None
    except Exception as e:
        print(f"Error al leer el archivo POSCAR: {e}")
        return None, None

    # --- 1. Parsear el POSCAR original ---
    try:
        # Línea 2: Factor de escala global (generalmente 1.0)
        global_scale_factor = float(lines[1].strip())

        # Líneas 3-5: Vectores de la celda unitaria
        primitive_lattice_vectors_raw = []
        for i in range(2, 5):
            vec = list(map(float, lines[i].strip().split()))
            primitive_lattice_vectors_raw.append(vec)
        primitive_lattice_matrix = np.array(primitive_lattice_vectors_raw)

        # Vectores de la celda unitaria reales (en Angstroms)
        primitive_lattice_matrix_angstrom = primitive_lattice_matrix * global_scale_factor

    except (IndexError, ValueError) as e:
        print(f"Error al parsear el POSCAR. Formato inesperado en las líneas de encabezado: {e}")
        return None, None

    # --- 2. Calcular los vectores de la supercelda ---
    supercell_lattice_matrix_angstrom = np.copy(primitive_lattice_matrix_angstrom)
    for i in range(3):
        supercell_lattice_matrix_angstrom[i, :] *= supercell_dim[i]

    # --- 3. Calcular el factor de escala de ALAMODE en Bohr ---
    # ALAMODE espera la longitud del primer vector de la supercelda en Bohr
    # para el factor de escala global de la red.
    alamode_scale_factor_angstrom = np.linalg.norm(supercell_lattice_matrix_angstrom[0])
    alamode_scale_factor_bohr = alamode_scale_factor_angstrom * ANGSTROM_TO_BOHR

    # --- 4. Calcular los vectores de la red normalizados para ALAMODE ---
    # Estos vectores son los vectores de la supercelda divididos por el
    # factor de escala en Bohr que acabamos de calcular.
    normalized_lattice_vectors = supercell_lattice_matrix_angstrom / alamode_scale_factor_bohr

    return alamode_scale_factor_bohr, normalized_lattice_vectors

# --- Bloque principal de ejecución ---
if __name__ == "__main__":
    print("--- Cálculo de Parámetros de Celda para ALAMODE ---")

    # Solicita la ruta del archivo POSCAR al usuario
    poscar_input_file = input("Introduce la ruta al archivo POSCAR de la celda unitaria (ej. POSCAR_unit): ")

    # Define las dimensiones de la supercelda
    # Modifica esto si necesitas una supercelda diferente a 2x2x2
    supercell_dimensions = (2, 2, 2)

    alamode_scale, normalized_vectors = get_alamode_cell_parameters(poscar_input_file, supercell_dimensions)

    if alamode_scale is not None and normalized_vectors is not None:
        print(f"\nResultados para la supercelda {supercell_dimensions[0]}x{supercell_dimensions[1]}x{supercell_dimensions[2]}:")
        print("\nPara tu archivo 'alm_suggest.in', deberías usar:")
        print("--------------------------------------------------------------------------------------------------")
        print(f"{alamode_scale:.16f}                     ! Factor de escala de la red en Bohr (Constante de red de la supercelda)")
        for vec in normalized_vectors:
            print(f"{vec[0]:.16f} {vec[1]:.16f} {vec[2]:.16f}")
        print("/&cutoff")
        print("!--------------------------------------------------------------------------------------------------")
        print("\nCopia y pega estas líneas en tu sección '&cell' de alm_suggest.in.")
    else:
        print("\nNo se pudieron calcular los parámetros de ALAMODE. Por favor, revisa el archivo POSCAR y el log de errores.")

    print("\nConsejo: Asegúrate de que las coordenadas atómicas en alm_suggest.in se generen ")
    print("en el bloque '/&position' para la SUPERCELDA y sean FRACCIONARIAS.")
    print("El script para convertir POSCAR direct a Cartesian (si es necesario) te dará las cartesianas.")
    print("Si necesitas las posiciones para alm_suggest.in, el script que te di antes")
    print("para construir el archivo alm_suggest.in a partir de SPOSCAR (cartesianas) y convertirlas a fraccionarias")
    print("ya tiene en cuenta la supercelda y las transformaciones de coordenadas.")
