#!/usr/bin/env python
import spglib
import numpy as np

# --- 1. Definir la celda unitaria para spglib ---
# Lattice vectors (Angstroms) for conventional FCC cell of Silicon
# Parametro de red 'a' para Si
a = 5.430

# Vectores de la red directa (celda convencional FCC)
lattice = np.array([
    [a, 0.0, 0.0],
    [0.0, a, 0.0],
    [0.0, 0.0, a]
])

# Posiciones de los átomos en coordenadas fraccionarias
# Para el Silicio (estructura diamante, 2 átomos por celda unitaria convencional)
positions = np.array([
    [0.0, 0.0, 0.0],
    [0.25, 0.25, 0.25]
])

# Tipos de átomos (pueden ser números, spglib no necesita los símbolos reales)
types = [0, 0] # 0 para Si (o cualquier número entero que represente el tipo)

# Ensamblar la celda en el formato spglib
cell = (lattice, positions, types)

# --- 2. Definir la ruta de alta simetría deseada (ej. FCC standard path) ---
# Puedes consultar tablas (como la de Setyawan y Curtarolo) para estas rutas.
# 'G' = Gamma, 'X', 'W', 'K', 'L' son puntos estándar para FCC.
kpath_string = "G-X-W-K-G|G-L" # Ejemplo de una ruta común para FCC

# --- 3. Obtener los puntos de alta simetría y la ruta ---
# num_points_per_segment es el numero de puntos que tendra cada segmento de la ruta.
# Si en anphon.in quieres 51 puntos, aquí pondrías 51.
num_points_per_segment = 51
path_data = spglib.get_kpoints_path(cell, kpath_string, num_points_per_segment)

# --- 4. Extraer los resultados ---
kpoints = path_data['points']
labels = path_data['labels']
path_segments = path_data['path']

print("--- Puntos de Alta Simetría Generados por spglib ---")
print("Total de puntos generados:", len(kpoints))

# Mostrar los puntos y sus etiquetas
for i, (point, label) in enumerate(zip(kpoints, labels)):
    if label: # Solo imprime la etiqueta si no está vacía
        print(f"  Punto: {label:<5} Coordenadas fraccionarias: [{point[0]:.6f}, {point[1]:.6f}, {point[2]:.6f}]")
    # else: # Para los puntos intermedios sin etiqueta
    #     print(f"  Punto {i+1:<3}: Coordenadas fraccionarias: [{point[0]:.6f}, {point[1]:.6f}, {point[2]:.6f}]")

print("\n--- Información de los Segmentos de la Ruta ---")
# El 'path_segments' te da los indices de inicio y fin de cada segmento en el array 'kpoints'
for i, segment in enumerate(path_segments):
    start_index, end_index = segment
    start_label = labels[start_index] if labels[start_index] else f"punto {start_index}"
    end_label = labels[end_index] if labels[end_index] else f"punto {end_index}"
    print(f"Segmento {i+1}: Del {start_label} (índice {start_index}) al {end_label} (índice {end_index}).")

# --- 5. Cómo usar esto en anphon.in ---
print("\n--- Formato para el bloque &kpoint en anphon.in (KPOINT_MODE = BAND) ---")
print("&kpoint")
print("  KPOINT_MODE = BAND")
# Puedes reconstruir los segmentos a partir de path_data['path'] y path_data['points']
# o simplemente tomar los puntos etiquetados directamente.

# Para el formato G 0.0 0.0 0.0 X 0.5 0.5 0.0 51
# Necesitas iterar sobre los segmentos y usar las coordenadas de los puntos etiquetados.
current_path = []
for i, (point, label) in enumerate(zip(kpoints, labels)):
    if label:
        current_path.append((label, point))
        if len(current_path) == 2:
            start_label, start_coords = current_path[0]
            end_label, end_coords = current_path[1]

            # Contar los puntos entre el inicio y el fin del segmento
            # Esto es un poco más complejo si path_data['path'] no da los indices exactos del segmento.
            # Una forma simple es asumir que cada segmento tiene num_points_per_segment.

            print(f"  {start_label} {start_coords[0]:.6f} {start_coords[1]:.6f} {start_coords[2]:.6f} "
                  f"{end_label} {end_coords[0]:.6f} {end_coords[1]:.6f} {end_coords[2]:.6f} {num_points_per_segment}")

            current_path = [current_path[1]] # Reiniciar para el siguiente segmento

# Si la ruta termina con un punto de simetría (ej. G-X-G), el ultimo G sera el inicio del siguiente o fin del ultimo.
# Este bucle simple asume pares. Ajustar para rutas continuas:
last_point_coords = None
last_point_label = None
for i, (point_coords, label) in enumerate(zip(kpoints, labels)):
    if label: # Es un punto de alta simetría
        if last_point_coords is not None:
            # Imprime el segmento anterior
            print(f"  {last_point_label} {last_point_coords[0]:.6f} {last_point_coords[1]:.6f} {last_point_coords[2]:.6f} "
                  f"{label} {point_coords[0]:.6f} {point_coords[1]:.6f} {point_coords[2]:.6f} {num_points_per_segment}")
        last_point_coords = point_coords
        last_point_label = label
    elif i > 0 and not label and labels[i-1]: # Punto intermedio del primer segmento
        # Esto es solo para capturar el primer segmento si no termina en un punto de etiqueta
        pass
# Este enfoque es más simple:
# spglib ya te da los puntos de alta simetria con etiquetas, el "num_points" es para los puntos intermedios.
# ALAMODE necesita el punto INICIAL, el punto FINAL y el NUMERO DE PUNTOS TOTAL en el segmento.
# La salida de spglib.get_kpoints_path ya te da TODOS los puntos (incluyendo los intermedios).

# Para el formato de ALAMODE:
# Simplemente iteramos sobre los segmentos proporcionados por spglib en 'path_data['path']'.
# Cada segmento [start_idx, end_idx] se refiere a los índices en 'path_data['points']'.
# Necesitamos la etiqueta y las coordenadas del punto de inicio y fin de cada segmento.
print("\n! Versión más robusta para el bloque &kpoint de anphon.in")
print("&kpoint")
print("  KPOINT_MODE = BAND")
for i, segment in enumerate(path_segments):
    start_idx, end_idx = segment
    start_point = kpoints[start_idx]
    end_point = kpoints[end_idx]

    # spglib puede tener etiquetas vacías para puntos intermedios.
    # Usamos la etiqueta del punto de inicio/fin si está disponible, sino una generica.
    start_label = labels[start_idx] if labels[start_idx] else f"P{start_idx}"
    end_label = labels[end_idx] if labels[end_idx] else f"P{end_idx}"

    # El número de puntos en el segmento es (end_idx - start_idx + 1)
    # Sin embargo, el 'num_points_per_segment' dado a spglib es el que ALAMODE espera por segmento.
    num_points_in_alamode_segment = num_points_per_segment

    print(f"  {start_label} {start_point[0]:.6f} {start_point[1]:.6f} {start_point[2]:.6f} "
          f"{end_label} {end_point[0]:.6f} {end_point[1]:.6f} {end_point[2]:.6f} "
          f"{num_points_in_alamode_segment}")
print("/")
