import numpy as np

def build_alm_input_file(sposcar_filepath: str, alm_header: str, alm_output_filepath: str):
    """
    Construye un archivo de entrada de ALAMODE a partir de un encabezado fijo
    y las posiciones de los átomos leídas de un SPOSCAR.

    Args:
        sposcar_filepath (str): Ruta al archivo SPOSCAR (generado por ASE, con coordenadas Cartesianas).
        alm_header (str): El contenido fijo del encabezado del archivo ALAMODE,
                          hasta justo antes del bloque '&position'.
        alm_output_filepath (str): Ruta donde se guardará el archivo ALAMODE final.
    """
    try:
        with open(sposcar_filepath, 'r') as f:
            sposcar_content = f.read()
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo SPOSCAR en {sposcar_filepath}")
        return

    sposcar_lines = sposcar_content.strip().split('\n')

    # --- 1. Parsear Vectores de Red del SPOSCAR (líneas 3 a 5, índice 2 a 4) ---
    lattice_vectors_str = [line.strip() for line in sposcar_lines[2:5]]
    lattice_matrix = []
    for vec_str in lattice_vectors_str:
        lattice_matrix.append(list(map(float, vec_str.split())))
    lattice_matrix = np.array(lattice_matrix) # Matriz cuyas filas son los vectores base

    # --- 2. Extraer Posiciones Atómicas Cartesianas ---
    try:
        cartesian_line_idx = sposcar_lines.index('Cartesian')
    except ValueError:
        print("Error: No se encontró 'Cartesian' en el SPOSCAR. Asegúrate de que las coordenadas sean cartesianas.")
        return

    # Las líneas de datos atómicos empiezan después de la línea 'Cartesian'
    atom_data_lines = [line.strip() for line in sposcar_lines[cartesian_line_idx + 1:] if line.strip()]

    # --- 3. Convertir Cartesianas a Fraccionarias y formatear ---
    formatted_positions = []
    inv_lattice_matrix = np.linalg.inv(lattice_matrix) # Inversa de la matriz base

    for line in atom_data_lines:
        coords_str = line.split()
        if len(coords_str) == 3:
            x, y, z = map(float, coords_str)

            # Convertir Cartesianas (x,y,z) a fraccionarias (u,v,w)
            # r_frac = r_cart . L_inv (donde r_cart es un vector fila, L_inv es la inversa de la matriz base)
            frac_coords = np.dot(np.array([x, y, z]), inv_lattice_matrix)

            # Formato para ALAMODE: "1 X Y Z" (1 para Si, tipo 1)
            formatted_positions.append(
                f"1 {frac_coords[0]:.16f} {frac_coords[1]:.16f} {frac_coords[2]:.16f}"
            )
        else:
            print(f"Advertencia: Saltando línea de posición mal formada en SPOSCAR: '{line}'")

    # --- 4. Construir el archivo ALAMODE final ---
    # Unir el encabezado, las posiciones formateadas y el delimitador final
    final_alm_content = alm_header.strip() + "\n" + "\n".join(formatted_positions) + "\n/"

    # Guardar el contenido en el archivo de salida
    with open(alm_output_filepath, "w") as f:
        f.write(final_alm_content)

    print(f"Archivo '{alm_output_filepath}' creado exitosamente con las coordenadas de la supercelda inyectadas y convertidas a fraccionarias.")


# --- Define el encabezado fijo para alm_suggest.in ---
# Este es el contenido de tu template hasta justo antes del bloque '&position'
ALM_SUGGEST_HEADER = """# ! alm_suggest.in (Formato adaptado a tu manual)
&general
  PREFIX = Si_displacements_patterns
  MODE = suggest
  NAT = 64
  NKD = 1
  KD = Si
  /

  &interaction

  NORDER = 2  # 1: armónico, 2: cúbico /&cell ! Orden máximo de las FC a considerar
# !--------------------------------------------------------------------------------------------------
20.5212                    # ! Factor de escala de la red en Bohr (Constante de red de la supercelda 2x2x2)
1.0 0.0 0.0                # ! Vectores de la red normalizados
0.0 1.0 0.0
0.0 0.0 1.0
/
&cutoff
# !--------------------------------------------------------------------------------------------------
Si-Si 5.0 3.0                 # ! Distancia de corte para interacciones armónicas (FC2) en Angstroms
#Si-Si-Si 3.0               # ! Distancia de corte para interacciones cúbicas (FC3) en Angstroms el formato no es así
/
&position
"""

# --- Ejecuta la función para construir alm_suggest.in ---
# Asegúrate de que tu archivo SPOSCAR (generado por el script de ASE)
# esté en el mismo directorio donde ejecutas este script.
build_alm_input_file(
    sposcar_filepath="SPOSCAR", # Asegúrate de que el archivo SPOSCAR exista aquí
    alm_header=ALM_SUGGEST_HEADER,
    alm_output_filepath="alm_suggest.in" # Sobrescribe o crea este archivo
)

# --- Ejemplos para alm_optimize.in y anphon.in (tendrías que adaptar sus headers) ---

# ALM_OPTIMIZE_HEADER = """# ! alm_optimize.in
# &general
#   PREFIX = Si_optimize
#   MODE = optimize
#   NAT = 64
#   NKD = 1
#   KD = Si
# /
# &interaction
#   NORDER = 2
# /
# &cell
# # !--------------------------------------------------------------------------------------------------
# 20.5212
# 1.0 0.0 0.0
# 0.0 1.0 0.0
# 0.0 0.0 1.0 /&cutoff
# # !--------------------------------------------------------------------------------------------------
# Si-Si 5.0 3.0
# #Si-Si-Si 3.0
# /
# &position
# """
# build_alm_input_file(
#     sposcar_filepath="SPOSCAR",
#     alm_header=ALM_OPTIMIZE_HEADER,
#     alm_output_filepath="alm_optimize.in"
# )

# ALM_ANPHON_HEADER = """! anphon.in
# &general
#   PREFIX = Si_LAMMPS_Phonons
#   MODE = phonon
#   NAT = 64
#   NKD = 1
#   KD = Si /&interaction
#   NORDER = 2
# /
# &cell
# # !--------------------------------------------------------------------------------------------------
# 20.5212
# 1.0 0.0 0.0
# 0.0 1.0 0.0
# 0.0 0.0 1.0
# /
# &cutoff
# # !--------------------------------------------------------------------------------------------------
# Si-Si 5.0 3.0
# #Si-Si-Si 3.0
# /
# &position
# """
# build_alm_input_file(
#     sposcar_filepath="SPOSCAR",
#     alm_header=ALM_ANPHON_HEADER,
#     alm_output_filepath="anphon.in"
# )

# Nota: Para 'alm_optimize.in' y 'anphon.in', aún necesitarás añadir manualmente
# los bloques finales con los parámetros específicos de optimización o fonones
# (como '&optimize' o '&phonon' y '&kpoint'). Este script solo se encarga
# de la parte de la estructura.
