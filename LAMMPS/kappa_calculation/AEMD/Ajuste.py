#!/usr/local/bin/python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
VOL = (5.431 * 50) * (2 * 5.431) * (2 * 5.431) * 1e-30 # m^3

# ==========================================
# 2. DEFINICIÓN DE FUNCIONES MATEMÁTICAS
# ==========================================
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
