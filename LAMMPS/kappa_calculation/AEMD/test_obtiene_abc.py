#!/usr/local/bin/python
"""
Este script es un ejemplo de la extracción de variables de un fichero de datos en
formato LAMMPS. 
"""
import numpy as np
import re

def get_lammps_vectors(file_path):
    """
    Lee un archivo de datos de LAMMPS y extrae los vectores de la red a, b, c.
    Maneja tanto celdas ortogonales como triclínicas.
    """
    # Inicialización de variables de borde
    data = {
        'xlo': 0.0, 'xhi': 0.0,
        'ylo': 0.0, 'yhi': 0.0,
        'zlo': 0.0, 'zhi': 0.0,
        'xy': 0.0, 'xz': 0.0, 'yz': 0.0
    }

    # Patrones para encontrar los límites y los factores de inclinación
    patterns = {
        'x': re.compile(r"^\s*(?P<lo>-?\d+\.?\d*)\s+(?P<hi>-?\d+\.?\d*)\s+xlo\s+xhi"),
        'y': re.compile(r"^\s*(?P<lo>-?\d+\.?\d*)\s+(?P<hi>-?\d+\.?\d*)\s+ylo\s+yhi"),
        'z': re.compile(r"^\s*(?P<lo>-?\d+\.?\d*)\s+(?P<hi>-?\d+\.?\d*)\s+zlo\s+zhi"),
        'tilt': re.compile(r"^\s*(?P<xy>-?\d+\.?\d*)\s+(?P<xz>-?\d+\.?\d*)\s+(?P<yz>-?\d+\.?\d*)\s+xy\s+xz\s+yz")
    }

    try:
        with open(file_path, 'r') as f:
            for line in f:
                # Buscar límites X, Y, Z
                for dim in ['x', 'y', 'z']:
                    match = patterns[dim].match(line)
                    if match:
                        data[f'{dim}lo'] = float(match.group('lo'))
                        data[f'{dim}hi'] = float(match.group('hi'))

                # Buscar factores de inclinación (tilt factors)
                tilt_match = patterns['tilt'].match(line)
                if tilt_match:
                    data['xy'] = float(tilt_match.group('xy'))
                    data['xz'] = float(tilt_match.group('xz'))
                    data['yz'] = float(tilt_match.group('yz'))

                # Si llegamos a Atoms o Masses, ya deberíamos tener la cabecera leída
                if "Atoms" in line or "Masses" in line:
                    break

        # Cálculo de las longitudes efectivas
        lx = data['xhi'] - data['xlo']
        ly = data['yhi'] - data['ylo']
        lz = data['zhi'] - data['zlo']

        # Construcción de los vectores según la convención de LAMMPS:
        # a = (xhi-xlo, 0, 0)
        # b = (xy, yhi-ylo, 0)
        # c = (xz, yz, zhi-zlo)

        a = np.array([lx, 0.0, 0.0])
        b = np.array([data['xy'], ly, 0.0])
        c = np.array([data['xz'], data['yz'], lz])

        return a, b, c

    except FileNotFoundError:
        print(f"Error: El archivo {file_path} no existe.")
        return None, None, None

# Ejemplo de uso con los datos proporcionados:
if __name__ == "__main__":
    # Creamos un archivo temporal para la prueba
    test_filename = "test_lammps.dat"
    with open(test_filename, "w") as f:
        f.write("# Cu3BHT Stacking AA\n\n")
        f.write("15 atoms\n3 atom types\n")
        f.write("0.000000000000 8.658944640366 xlo xhi\n")
        f.write("0.000000000000 7.522840153334 ylo yhi\n")
        f.write("0.000000000000 3.183817710818 zlo zhi\n")
        f.write("4.324263184156 -0.016658943914 1.195978650157 xy xz yz\n")
        f.write("Masses\n1 32.06\n")

    va, vb, vc = get_lammps_vectors(test_filename)

    print("Vector a:", va)
    print("Vector b:", vb)
    print("Vector c:", vc)
