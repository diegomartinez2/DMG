import numpy as np
from scipy.integrate import simps # o trapz para integración
import os

# =======================================================
# 1. Parámetros de la Simulación y Constantes Físicas
# =======================================================

# Parámetros usados en LAMMPS (debe coincidir con el input file)
T = 80.0          # Temperatura del sistema (K)
V = 1475.0        # Volumen de la caja de simulación (A^3) - Ejemplo, usar el valor real de su simulación
dt_lammps = 0.001 # Paso de tiempo de la simulación (ps)
#dt_lammps = 0.0005

# Constantes Físicas (en unidades de LAMMPS: eV, A, ps)
# k_B (Constante de Boltzmann) en eV/K
kB_eV_K = 8.617333262145e-5

# Factor de Conversión de Unidades: (eV^2 / (A^4 ps K)) a W/(m K)
# Unidades LAMMPS: (Flujo)^2 * tiempo = (eV/A^2/ps)^2 * ps * A^3 * ps = eV^2 / (A ps)
# Factor = (eV^2 / (A^4 ps K)) * (A^3 ps) / (eV^2 / (A ps))
# Factor de conversión: 1 eV^2 / (A^4 ps) = 1.0e-29 W / (m K)
# Más fácil: El resultado de la fórmula de GK en unidades "metal" es eV^2 / (A ps K)
# Factor de conversión (eV^2 / (A ps K)) a W/(m K):
FACTOR_CONVERSION = 1.5975e+10

# Tiempo de truncamiento para la integral (en ps).
# Escoger donde la ACF se acerca a cero o se vuelve ruidosa.
tau_max = 50.0
#tau_max = 100.0
# =======================================================
# 2. Lectura y Ajuste de Datos
# =======================================================

input_file = 'J0Jt.dat'

# Verificar si el archivo existe
if not os.path.exists(input_file):
    print(f"Error: Archivo {input_file} no encontrado.")
    exit()

# Leer los datos, saltando las líneas de encabezado de LAMMPS
# Se asume que las primeras 4 líneas son comentarios/encabezado
data = np.loadtxt(input_file, skiprows=4)

# Columna 0: Step, Columna 1-3: ACF (JxJx, JyJy, JzJz)
# El tiempo (t) se calcula a partir del número de paso y el paso de tiempo
time_steps = data[:, 0]
time = time_steps * dt_lammps # Tiempo en ps

# ACF del flujo de calor, normalizada por el número de pasos promediados (LAMMPS ya lo hace)
acf_x = data[:, 1]
acf_y = data[:, 2]
acf_z = data[:, 3]

# # Cargar múltiples archivos J0Jt.dat
# input_files = ['J0Jt_12345.dat', 'J0Jt_23456.dat', 'J0Jt_34567.dat']
# acf_x_avg = np.zeros_like(np.loadtxt(input_files[0], skiprows=4)[:, 1])
# for f in input_files:
#     data = np.loadtxt(f, skiprows=4)
#     acf_x_avg += data[:, 1]
# acf_x_avg /= len(input_files)
# # Procede con integración usando acf_x_avg

# Recortar datos hasta el tiempo de truncamiento (tau_max)
idx_max = np.where(time >= tau_max)[0]
if len(idx_max) > 0:
    idx_max = idx_max[0]
else:
    idx_max = len(time)

time_trunc = time[:idx_max]
acf_x_trunc = acf_x[:idx_max]
acf_y_trunc = acf_y[:idx_max]
acf_z_trunc = acf_z[:idx_max]

# =======================================================
# 3. Integración de la Función de Autocorrelación
# =======================================================

# Usar la regla de Simpson (simps) o Trapecio (trapz) para la integración numérica
# El resultado de la integración es \int_0^\tau <J(t)J(0)> dt
Int_x = simps(acf_x_trunc, time_trunc)
Int_y = simps(acf_y_trunc, time_trunc)
Int_z = simps(acf_z_trunc, time_trunc)

# Promedio de las tres direcciones
Int_avg = (Int_x + Int_y + Int_z) / 3.0

# =======================================================
# 4. Cálculo de la Conductividad Térmica (kappa)
# =======================================================

# Fórmula de Green-Kubo: kappa = (1 / (3 V k_B T^2)) * \int_0^\tau <J(t)J(0)> dt
# El flujo de calor en LAMMPS (compute heat/flux) ya incluye V, así que la fórmula es:
# kappa = (1 / (3 k_B T^2)) * \int_0^\tau <J(t)J(0)> dt
# Donde <J(t)J(0)> es el resultado de la correlación de compute heat/flux
# El output de LAMMPS (J0Jt.dat) tiene unidades de eV^2/(A^2 ps^2) * A^6 = eV^2 A^4 / ps^2

# Calculando el factor (1 / (3 * k_B * T^2))
#GK_prefactor = 1.0 / (3.0 * kB_eV_K * T**2) #LAMMPS ya hace la salida de los Jx Jy Jz por volumen (pero hay que revisar esto)
GK_prefactor = V / (3.0 * kB_eV_K * T**2)

# Conductividad térmica en las unidades de LAMMPS (eV^2 / (A ps K))
kappa_x_LAMMPS = GK_prefactor * Int_x
kappa_y_LAMMPS = GK_prefactor * Int_y
kappa_z_LAMMPS = GK_prefactor * Int_z
kappa_avg_LAMMPS = GK_prefactor * Int_avg

# 5. Conversión a Unidades Internacionales (SI): W/(m K)
kappa_x_SI = kappa_x_LAMMPS * FACTOR_CONVERSION
kappa_y_SI = kappa_y_LAMMPS * FACTOR_CONVERSION
kappa_z_SI = kappa_z_LAMMPS * FACTOR_CONVERSION
kappa_avg_SI = kappa_avg_LAMMPS * FACTOR_CONVERSION

# =======================================================
# 6. Resultados
# =======================================================

print(f"\n--- Resultados del Cálculo de Green-Kubo ---")
print(f"Temperatura (T): {T:.2f} K")
print(f"Volumen de Simulación (V): {V:.2f} Å^3")
print(f"Tiempo de Truncamiento (tau_max): {tau_max:.1f} ps")
print("-" * 40)
print(f"Integral de ACF (X): {Int_x:.2e} [Unidades LAMMPS]")
print(f"Integral de ACF (Y): {Int_y:.2e} [Unidades LAMMPS]")
print(f"Integral de ACF (Z): {Int_z:.2e} [Unidades LAMMPS]")
print("-" * 40)
print(f"Conductividad Térmica (kappa) en W/(m K):")
print(f"kappa_x: {kappa_x_SI:.3f}")
print(f"kappa_y: {kappa_y_SI:.3f}")
print(f"kappa_z: {kappa_z_SI:.3f}")
print(f"kappa_promedio: {kappa_avg_SI:.3f} W/(m K) 🌡️")

# Opcional: Graficar la ACF y su Integral
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(time_trunc, acf_x_trunc, label=r'$\langle J_x(t) J_x(0) \rangle$')
plt.plot(time_trunc, acf_y_trunc, label=r'$\langle J_y(t) J_y(0) \rangle$')
plt.plot(time_trunc, acf_z_trunc, label=r'$\langle J_z(t) J_z(0) \rangle$')
plt.axhline(0, color='grey', linestyle='--')
plt.title('Función de Autocorrelación del Flujo de Calor')
plt.xlabel('Tiempo (ps)')
plt.ylabel('ACF del Flujo de Calor [Unidades LAMMPS]')
plt.legend()
plt.grid(True)
plt.show()

# **Cambio: Grafica la integral acumulativa para verificar convergencia (plateau)**
cum_int_x = np.cumsum(acf_x_trunc) * dt_lammps  # Integral acumulativa
cum_int_y = np.cumsum(acf_y_trunc) * dt_lammps  # Integral acumulativa
cum_int_z = np.cumsum(acf_z_trunc) * dt_lammps  # Integral acumulativa
plt.figure(figsize=(10, 6))
plt.plot(time_trunc, cum_int_x, label='Integral Acumulativa X')
plt.plot(time_trunc, cum_int_y, label='Integral Acumulativa Y')
plt.plot(time_trunc, cum_int_z, label='Integral Acumulativa Z')
plt.title('Integral Acumulativa de ACF para Convergencia')
plt.xlabel('Tiempo (ps)')
plt.ylabel('Integral [Unidades LAMMPS]')
plt.legend()
plt.grid(True)
plt.show()
