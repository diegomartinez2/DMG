import numpy as np
from scipy.integrate import simps
import os

def read_lammps_params(filename="parametros.txt"):
    """Lee T, V y dt del archivo parametros.txt generado por LAMMPS."""
    parametros = {}
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Error: Archivo de parámetros {filename} no encontrado.")

    with open(filename, 'r') as f:
        for line in f:
            if ':' in line:
                key, value = line.split(':')
                # Limpiar la cadena y convertir a float
                parametros[key.strip()] = float(value.strip())

    return parametros

def calculate_kappa(parametros, input_file='J0Jt.dat'):
    """Calcula la conductividad térmica usando el método Green-Kubo."""

    T = parametros['T']
    V = parametros['V']
    dt_lammps = parametros['dt']
    tau_max = 50.0 # Tiempo de truncamiento (ps), puede ser fijo o pasado como arg.
#tau_max = 100.0

    # Constantes Físicas y Conversión
    kB_eV_K = 8.617333262145e-5 # Constante de Boltzmann en eV/K
    # Factor de Conversión de Unidades: (eV^2 / (A^4 ps K)) a W/(m K)
    # Unidades LAMMPS: (Flujo)^2 * tiempo = (eV/A^2/ps)^2 * ps * A^3 * ps = eV^2 / (A ps)
    # Factor = (eV^2 / (A^4 ps K)) * (A^3 ps) / (eV^2 / (A ps))
    # Factor de conversión: 1 eV^2 / (A^4 ps) = 1.0e-29 W / (m K)
    # Más fácil: El resultado de la fórmula de GK en unidades "metal" es eV^2 / (A ps K)
    # Factor de conversión (eV^2 / (A ps K)) a W/(m K):
    FACTOR_CONVERSION = 1.5975e+10
    #FACTOR_CONVERSION = 1.602176634e+10 # Factor de (eV*A^2 / (K*ps)) a W/(m K)

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
    acf_x_trunc = acf_x[:idx_max]
    acf_y_trunc = acf_y[:idx_max]
    acf_z_trunc = acf_z[:idx_max]

    # Integración de la ACF
    Int_x = simps(acf_x[:idx_max], time_trunc)
    Int_y = simps(acf_y[:idx_max], time_trunc)
    Int_z = simps(acf_z[:idx_max], time_trunc)
    # Promedio de las tres direcciones
    Int_avg = (Int_x + Int_y + Int_z) / 3.0
    # Fórmula de Green-Kubo (asumiendo compute heat/flux sin V en el prefactor)
    GK_prefactor = 1.0 / (3.0 * kB_eV_K * T**2)
    GK_prefactor = V / (3.0 * kB_eV_K * T**2)

    # Conductividad térmica en W/(m K)
    kappa_x_SI = GK_prefactor * Int_x * FACTOR_CONVERSION
    kappa_y_SI = GK_prefactor * Int_y * FACTOR_CONVERSION
    kappa_z_SI = GK_prefactor * Int_z * FACTOR_CONVERSION
    kappa_avg_SI = (kappa_x_SI + kappa_y_SI + kappa_z_SI) / 3.0
    kappa_avg2_SI = GK_prefactor * Int_avg * FACTOR_CONVERSION
    print(f"\n--- Resultados de Conductividad Térmica ---")
    print(f"Temperatura: {T:.2f} K | Volumen: {V:.2f} A^3")
    print(f"Tiempo de Correlación: {time_trunc[-1]:.1f} ps")
    print("-" * 40)
    print(f"Kappa_x: {kappa_x_SI:.3f} W/(m K)")
    print(f"Kappa_y: {kappa_y_SI:.3f} W/(m K)")
    print(f"Kappa_z: {kappa_z_SI:.3f} W/(m K)")
    print("-" * 40)
    print(f"Kappa Promedio: {kappa_avg_SI:.3f} W/(m K) 🌡️")
    #print(f"Kappa Promedio(2): {kappa_avg2_SI:.3f} W/(m K) 🌡️")
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
    print(f"kappa_promedio: {kappa_avg2_SI:.3f} W/(m K) 🌡️")

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


if __name__ == "__main__":
    try:
        sim_params = read_lammps_params()
        calculate_kappa(sim_params)
    except Exception as e:
        print(f"Fallo en el post-procesamiento: {e}")
