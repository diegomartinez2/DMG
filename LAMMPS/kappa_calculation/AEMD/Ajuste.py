#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Ajuste.py
#
#  Copyright 2026 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
#
"""
Este script ajusta los resultados del cálculo AEMD para obtener la conductividad
térmica. Puede leer los datos del volumen del archivo de LAMMPS de datos.
"""
# ---------------------------
# Importación de los módulos
# ---------------------------
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re

# ==========================================
# 1. PARÁMETROS CONFIGURABLES (Variables)
# ==========================================
L = 135.775       # Distancia del sistema
T_MAX = 400.0     # Temperatura bloque caliente
T_MIN = 200.0     # Temperatura bloque frío
N_TERMS = 20      # Número de términos en la serie de Fourier
X_OFFSET = 20000  # Desplazamiento en el eje X ($1 - 20000)

# Constantes para el cálculo final (kbarr)
KB = 1.38e-23     # Constante de Boltzmann
N_PARTICLES = 1600
DIM = 3
VOL = (50 * 5.431) * (2 * 5.431) * (2 * 5.431) * 1e-30 # m^3

# ==========================================
# 2. DEFINICIÓN DE FUNCIONES MATEMÁTICAS
# ==========================================
def Volumen_celda(arg):
    """a,b,c se obtienen de la celda triclinica que se define con
    xlo xhi, ylo yhi, zlo zhi y los factores de tilt xy xz yz """
    a = np.array([8.658944640366, 0.0, 0.0])                        #ejemplo
    b = np.array([4.324263184156, 7.522840153334, 0.0])             #ejemplo
    c = np.array([-0.016658943914, 1.195978650157, 3.183817710818]) #ejemplo

    vol = np.abs(np.dot(a, np.cross(b, c)))
    print(f"Volumen: {vol:.3f} Å³")
    print(f"Vol/átomo: {vol/15:.3f} Å³/átomo")
    pass

def get_lammps_vectors(file_path="log.lamps"):
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

def calcular_volumen_super(a, b, c, nx=1, ny=1, nz=1):
    """Versión generalizada para el calculo del volumen de una supercelda
    uso:
           v_unit, v_super = calcular_volumen_super(a, b, c, nx, ny, nz)"""
    vol_unit = np.abs(np.dot(a, np.cross(b, c)))
    vol_super = vol_unit * nx * ny * nz
    print(f"Volumen supercelda: {vol_super:.3f} Å³")
    return vol_unit, vol_super


def extraer_parametros_lammps(file_path="log.lammps"):
    """
    Analiza un archivo log de LAMMPS para extraer dimensiones, volumen,
    número de átomos y temperaturas del pulso AEMD.
    """
    data = {
        "volumen_ang3": None,
        "longitud_x": None,
        "longitud_y": None,
        "longitud_z": None,
        "num_atomos": None,
        "t_caliente": None,
        "t_fria": None
    }

    try:
        with open(file_path, 'r') as f:
            content = f.read()

        # 1. Extraer dimensiones de la caja (buscamos la última caja mencionada tras replicate)
        # Formato: triclinic box = (0 0 0) to (x y z)
        box_match = re.findall(r"triclinic box = \(0 0 0\) to \(([\d\.]+) ([\d\.]+) ([\d\.]+)\)", content)
        if box_match:
            # Tomamos la última coincidencia (la del sistema ya replicado)
            lx, ly, lz = map(float, box_match[-1])
            data["longitud_x"] = lx
            data["longitud_y"] = ly
            data["longitud_z"] = lz
            # Nota: Para cajas triclínicas, el volumen es lx * ly * lz independientemente del tilt
            data["volumen_ang3"] = lx * ly * lz

        # 2. Extraer número de átomos (después de la replicación)
        atoms_match = re.findall(r"(\d+) atoms", content)
        if atoms_match:
            # Tomamos el último valor reportado antes de las simulaciones
            data["num_atomos"] = int(atoms_match[-1])

        # 3. Extraer temperaturas (buscando las definiciones de variables en el script)
        t_hot_match = re.search(r"variable T_hot_pulse equal ([\d\.]+)", content)
        t_cold_match = re.search(r"variable T_cold_pulse equal ([\d\.]+)", content)

        if t_hot_match:
            data["t_caliente"] = float(t_hot_match.group(1))
        if t_cold_match:
            data["t_fria"] = float(t_cold_match.group(1))

        return data

    except FileNotFoundError:
        return "Error: El archivo no fue encontrado."
    except Exception as e:
        return f"Error al procesar el archivo: {e}"


    resultados = data

    if isinstance(resultados, dict):
        print("--- Parámetros Extraídos ---")
        print(f"Número de átomos: {resultados['num_atomos']}")
        print(f"Longitud X: {resultados['longitud_x']} Å")
        print(f"Volumen: {resultados['volumen_ang3']:.2f} Å³")
        print(f"Temperatura Caliente: {resultados['t_caliente']} K")
        print(f"Temperatura Fría: {resultados['t_fria']} K")
    else:
        print(resultados)

def get_constants(n):
    """Calcula los coeficientes A, B y k basados en n."""
    k_n = (np.pi * n) / L
    k1_n = np.pi * n

    # Coeficiente B(n) según el script original
    # B(n) = ((1/k1(n))*(-cos(k1(n))*(Tmin-Tmax)-Tmax+Tmin))
    b_n = (1.0 / k1_n) * (-np.cos(k1_n) * (T_MIN - T_MAX) - T_MAX + T_MIN)

    return k_n, k1_n, b_n

def f_n(x, n, k2):
    """Cálculo de un término individual de la serie."""
    kn, k1n, bn = get_constants(n)
    # Coeficiente C(n) integrado directamente en la función de decaimiento
    term_c = (-2 * bn + 2 * bn * np.cos(k1n)) * (1.0 / k1n)
    return np.exp(-k2 * x * (kn**2)) * term_c

def f_total(x, k2):
    """Sumatoria F0(x) de n=1 hasta N_TERMS."""
    total = 0
    for n in range(1, N_TERMS + 1):
        total += f_n(x, n, k2)
    return total

# ==========================================
# 3. CARGA DE DATOS Y AJUSTE (FIT)
# ==========================================
# Nota: Reemplaza 'out1' por la ruta de tu archivo.
# Si no existe, el script generará datos sintéticos para demostración.
try:
    data = np.loadtxt('out1')
    x_data = data[:, 0] - X_OFFSET
    y_data = data[:, 3] # Columna 4 en Gnuplot es índice 3 en Python
#    y_data = data[:, 1] # Columna 2 en Gnuplot es índice 1 en Python
except FileNotFoundError:
    print("Archivo 'out1' no encontrado. Generando datos de ejemplo...")
    x_data = np.linspace(0, 2000000, 100)
    y_data = f_total(x_data, 0.00002) + np.random.normal(0, 0.05, 100)

# Ajuste (Fit) de la variable K2
# p0 es el valor inicial (K2=0.00001 en tu script)
popt, pcov = curve_fit(f_total, x_data, y_data, p0=[0.00001])
k2_fit = popt[0]

# ==========================================
# 4. CÁLCULOS FINALES Y SALIDA
# ==========================================
# kbarr = (K2 * 1e-5 * 3 * 1600 * 1.38e-23) / Volumen
kbarr = (k2_fit * 1e-5 * DIM * N_PARTICLES * KB) / VOL

print("-" * 30)
print(f"RESULTADOS:")
print(f"K2 optimizado: {k2_fit:.6e}")
print(f"Conductividad Térmica (kbarr): {kbarr:.6f}")
print("-" * 30)

# ==========================================
# 5. GRÁFICA "BONITA"
# ==========================================
plt.figure(figsize=(10, 6))
plt.style.use('seaborn-v0_8-whitegrid') # Estilo limpio

plt.scatter(x_data, y_data, label='Datos (out1)', alpha=0.5, color='gray', s=10)
plt.plot(x_data, f_total(x_data, k2_fit), label=f'Ajuste F0(x)\n$K_2$ = {k2_fit:.2e}', color='firebrick', lw=2)

plt.title('Ajuste de Difusividad Térmica', fontsize=14)
plt.xlabel('Tiempo / Distancia desplazada ($x - 20000$)', fontsize=12)
plt.ylabel('Amplitud / Temperatura', fontsize=12)
plt.legend(frameon=True)
plt.tight_layout()

plt.show()
