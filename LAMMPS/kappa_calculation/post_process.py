import numpy as np
from scipy.integrate import simps
import os

def read_lammps_params(filename="params.txt"):
    """Lee T, V y dt del archivo params.txt generado por LAMMPS."""
    params = {}
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Error: Archivo de parámetros {filename} no encontrado.")

    with open(filename, 'r') as f:
        for line in f:
            if ':' in line:
                key, value = line.split(':')
                # Limpiar la cadena y convertir a float
                params[key.strip()] = float(value.strip())

    return params

def calculate_kappa(params, input_file='J0Jt.dat'):
    """Calcula la conductividad térmica usando el método Green-Kubo."""

    T = params['T']
    V = params['V']
    dt_lammps = params['dt']
    tau_max = 50.0 # Tiempo de truncamiento (ps), puede ser fijo o pasado como arg.

    # Constantes Físicas y Conversión
    kB_eV_K = 8.617333262145e-5 # Constante de Boltzmann en eV/K
    FACTOR_CONVERSION = 1.602176634e+10 # Factor de (eV*A^2 / (K*ps)) a W/(m K)

    if not os.path.exists(input_file):
        print(f"Error: Archivo de ACF {input_file} no encontrado.")
        return None

    # Leer los datos de la ACF
    data = np.loadtxt(input_file, skiprows=4)
    time_steps = data[:, 0]
    time = time_steps * dt_lammps # Tiempo en ps

    # ACF del flujo de calor (columnas 1, 2, 3)
    acf_x, acf_y, acf_z = data[:, 1], data[:, 2], data[:, 3]

    # Recortar datos
    idx_max = np.where(time >= tau_max)[0]
    if len(idx_max) > 0:
        idx_max = idx_max[0]
    else:
        idx_max = len(time)

    time_trunc = time[:idx_max]

    # Integración de la ACF
    Int_x = simps(acf_x[:idx_max], time_trunc)
    Int_y = simps(acf_y[:idx_max], time_trunc)
    Int_z = simps(acf_z[:idx_max], time_trunc)

    # Fórmula de Green-Kubo (asumiendo compute heat/flux sin V en el prefactor)
    GK_prefactor = 1.0 / (3.0 * kB_eV_K * T**2)

    # Conductividad térmica en W/(m K)
    kappa_x_SI = GK_prefactor * Int_x * FACTOR_CONVERSION
    kappa_y_SI = GK_prefactor * Int_y * FACTOR_CONVERSION
    kappa_z_SI = GK_prefactor * Int_z * FACTOR_CONVERSION
    kappa_avg_SI = (kappa_x_SI + kappa_y_SI + kappa_z_SI) / 3.0

    print(f"\n--- Resultados de Conductividad Térmica ---")
    print(f"Temperatura: {T:.2f} K | Volumen: {V:.2f} A^3")
    print(f"Tiempo de Correlación: {time_trunc[-1]:.1f} ps")
    print("-" * 40)
    print(f"Kappa_x: {kappa_x_SI:.3f} W/(m K)")
    print(f"Kappa_y: {kappa_y_SI:.3f} W/(m K)")
    print(f"Kappa_z: {kappa_z_SI:.3f} W/(m K)")
    print("-" * 40)
    print(f"Kappa Promedio: {kappa_avg_SI:.3f} W/(m K) 🌡️")

if __name__ == "__main__":
    try:
        sim_params = read_lammps_params()
        calculate_kappa(sim_params)
    except Exception as e:
        print(f"Fallo en el post-procesamiento: {e}")
