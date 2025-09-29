#!/usr/bin/env python
import numpy as np

def convert_poscar_direct_to_cartesian(input_filepath: str, output_filepath: str):
    """
    Transforma un archivo POSCAR de coordenadas 'Direct' (fraccionales)
    a un archivo POSCAR con coordenadas 'Cartesian' (absolutas).

    Args:
        input_filepath (str): Ruta al archivo POSCAR de entrada con coordenadas 'Direct'.
        output_filepath (str): Ruta donde se guardará el nuevo archivo POSCAR con coordenadas 'Cartesian'.
    """
    try:
        with open(input_filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo de entrada en {input_filepath}")
        return

    # --- 1. Leer las primeras 8 líneas (encabezado) ---
    # Línea 1: Comentario/Nombre del sistema
    # Línea 2: Factor de escala global
    # Líneas 3-5: Vectores de la celda de la red (a1, a2, a3)
    # Línea 6: Nombres de los tipos de átomos
    # Línea 7: Conteo de átomos por tipo
    # Línea 8: Tipo de coordenadas ('Direct' o 'Cartesian')

    header_lines = lines[:8]
    output_lines = header_lines[:] # Copiar las líneas del encabezado para la salida

    # Obtener el factor de escala global
    scale_factor = float(lines[1].strip())

    # Obtener los vectores de la red
    lattice_vectors = []
    for i in range(2, 5):
        vec = list(map(float, lines[i].strip().split()))
        lattice_vectors.append(vec)
    lattice_matrix = np.array(lattice_vectors) * scale_factor # Multiplicar por el factor de escala

    # Asegurarse de que el tipo de coordenadas es 'Direct'
    coords_type_line_idx = 7
    current_coords_type = lines[coords_type_line_idx].strip().lower()

    if current_coords_type not in ['direct', 'cartesian']:
        # Handle cases where Selective Dynamics or another line might be present
        # In VASP POSCAR, if Selective Dynamics is present, it's usually line 8,
        # and the coordinate type is line 9.
        # Let's try to find it dynamically by looking for 'Direct' or 'Cartesian'
        found_coords_line = False
        for i in range(7, 10): # Check a few lines just in case
            if i < len(lines):
                if lines[i].strip().lower() == 'direct':
                    coords_type_line_idx = i
                    current_coords_type = 'direct'
                    found_coords_line = True
                    break
                elif lines[i].strip().lower() == 'cartesian':
                    print("Advertencia: El archivo POSCAR de entrada ya está en formato 'Cartesian'. No se realizará la conversión.")
                    return
        if not found_coords_line:
             print("Error: No se pudo determinar el tipo de coordenadas (esperado 'Direct' o 'Cartesian').")
             print(f"Línea 8 (o cercana) encontrada: '{lines[7].strip()}'")
             return


    if current_coords_type == 'cartesian':
        print("Advertencia: El archivo POSCAR de entrada ya está en formato 'Cartesian'. No se realizará la conversión.")
        return

    # --- 2. Modificar la línea del tipo de coordenadas a 'Cartesian' ---
    output_lines[coords_type_line_idx] = "Cartesian\n"

    # --- 3. Leer las coordenadas atómicas y convertirlas ---
    num_atoms_per_type = list(map(int, lines[6].strip().split()))
    total_atoms = sum(num_atoms_per_type)

    # Las coordenadas de los átomos comienzan después de la línea del tipo de coordenadas
    start_coords_idx = coords_type_line_idx + 1

    # Verificar si hay una línea de 'Selective Dynamics'
    selective_dynamics_line = ""
    if lines[start_coords_idx].strip().lower() == 'selective dynamics':
        selective_dynamics_line = lines[start_coords_idx]
        output_lines.insert(coords_type_line_idx + 1, selective_dynamics_line) # Mantenerla si existe
        start_coords_idx += 1 # Ajustar el inicio de las coordenadas

    direct_coords_lines = lines[start_coords_idx : start_coords_idx + total_atoms]

    if len(direct_coords_lines) != total_atoms:
        print(f"Error: El número de líneas de coordenadas ({len(direct_coords_lines)}) no coincide con el número total de átomos ({total_atoms}).")
        return

    converted_coords = []
    for line in direct_coords_lines:
        parts = line.strip().split()
        direct_coords = np.array(list(map(float, parts[:3]))) # Solo las 3 primeras columnas

        # Convertir coordenadas directas (fraccionales) a cartesianas
        # r_cart = u*a1 + v*a2 + w*a3
        # Esto es equivalente a np.dot(direct_coords, lattice_matrix) si direct_coords es (u,v,w)
        cartesian_coords = np.dot(direct_coords, lattice_matrix)

        # Formatear la línea de salida, manteniendo cualquier información extra (Selective Dynamics, comentarios)
        extra_info = " ".join(parts[3:]) if len(parts) > 3 else ""
        converted_coords.append(
            f" {cartesian_coords[0]:.16f} {cartesian_coords[1]:.16f} {cartesian_coords[2]:.16f} {extra_info}".strip() + "\n"
        )

    # --- 4. Escribir el nuevo archivo POSCAR ---
    # Reemplazar las antiguas líneas de coordenadas con las nuevas
    # Primero, eliminar las antiguas líneas de coordenadas si Selective Dynamics no estaba presente
    # Si Selective Dynamics sí estaba, ya la insertamos y solo debemos añadir las coordenadas después.

    # En este enfoque, reconstruimos el archivo desde el encabezado y las nuevas coordenadas
    final_content = "".join(output_lines[:coords_type_line_idx + 1]) # Encabezado hasta 'Cartesian'
    if selective_dynamics_line:
        final_content += selective_dynamics_line # Añadir Selective Dynamics si estaba presente
    final_content += "".join(converted_coords) # Añadir las nuevas coordenadas

    with open(output_filepath, 'w') as f:
        f.write(final_content)

    print(f"Conversión completada: '{input_filepath}' (Direct) convertido a '{output_filepath}' (Cartesian).")

# --- Cómo usar el script ---
if __name__ == "__main__":
    # Nombre del archivo POSCAR de entrada (ej. tu_poscar_direct.vasp)
    input_poscar = "AA_relaxed.vasp"
    # Nombre del archivo POSCAR de salida
    output_poscar = "AA_relaxed.POSCAR_cartesian"

    # Ejecutar la función de conversión
    convert_poscar_direct_to_cartesian(input_poscar, output_poscar)

    print("\nPara usar este script:")
    print("1. Guarda el código en un archivo llamado 'poscar_direct_to_cartesian.py'.")
    print("2. Asegúrate de que tu archivo POSCAR de entrada (ej. 'POSCAR_direct') esté en el mismo directorio.")
    print("3. Ejecuta el script desde la terminal: 'python poscar_direct_to_cartesian.py'")
    print("   Se generará un nuevo archivo llamado 'POSCAR_cartesian'.")
    print("\nNota: El script intenta manejar la línea 'Selective Dynamics' si está presente.")
