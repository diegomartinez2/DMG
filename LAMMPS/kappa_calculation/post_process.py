#!/usr/local/bin/python
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
    GK_prefactor = 1.0 / (3.0 * kB_eV_K * T**2) #para los datos de 'J0Jt.dat'
    #GK_prefactor = V / (3.0 * kB_eV_K * T**2) #para los datos guardados en 'J0Jt_12345.dat' como Js (sin ser Js/vol)

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

def calculate_kappa_promedio(parametros, input_files=['J0Jt_12345.dat']):
    """Calcula la conductividad térmica usando el método Green-Kubo.

    Args:
        parametros (dict): Diccionario con T (K), V (Å^3), dt (ps), tau_max (ps).
        input_files (list): Lista de archivos J0Jt_*.dat para promediar.

    Returns:
        dict: Resultados de kappa (W/(m K)) y ACF integrales.
    """
    # Extraer parámetros
    T = parametros.get('T', 80.0)
    V = parametros.get('V', None)
    dt_lammps = parametros.get('dt', 0.001)
    tau_max = parametros.get('tau_max', 100.0)  # Default para materiales complejos

    # Validar parámetros
    if not all(isinstance(x, (int, float)) and x > 0 for x in [T, V, dt_lammps, tau_max]):
        print("Error: Parámetros T, V, dt, tau_max deben ser positivos.")
        return None

    # Constantes físicas
    kB_eV_K = 8.617333262145e-5  # Boltzmann en eV/K
    FACTOR_CONVERSION = 1.60217662e13 / 1e10  # eV^2/(Å·ps·K) a W/(m·K)
    # Derivación: 1 eV = 1.60217662e-19 J, 1 Å = 1e-10 m, 1 ps = 1e-12 s
    # ACF units: (eV/Å^2/ps)^2·Å^6·ps = eV^2·Å^4/ps
    # kappa = V * Int / (3 kB T^2) -> eV^2·Å^4/(ps·K) * Å^3 / K = eV^2·Å^7/(ps·K^2)
    # Convert to W/(m·K): factor = (1.60217662e-19)^2 / (1e-10)^7 / (1e-12) = 1.60217662e13

    # Verificar archivos
    for f in input_files:
        if not os.path.exists(f):
            print(f"Error: Archivo {f} no encontrado.")
            return None

    # Leer y promediar ACFs
    acf_x_avg, acf_y_avg, acf_z_avg = None, None, None
    for i, f in enumerate(input_files):
        data = np.loadtxt(f, skiprows=4)
        if i == 0:
            time_steps = data[:, 0]
            acf_x_avg = data[:, 1]
            acf_y_avg = data[:, 2]
            acf_z_avg = data[:, 3]
        else:
            acf_x_avg += data[:, 1]
            acf_y_avg += data[:, 2]
            acf_z_avg += data[:, 3]
    acf_x_avg /= len(input_files)
    acf_y_avg /= len(input_files)
    acf_z_avg /= len(input_files)
    time = time_steps * dt_lammps

    # Recortar datos hasta tau_max
    idx_max = np.where(time >= tau_max)[0]
    idx_max = idx_max[0] if len(idx_max) > 0 else len(time)
    time_trunc = time[:idx_max]
    acf_x_trunc = acf_x_avg[:idx_max]
    acf_y_trunc = acf_y_avg[:idx_max]
    acf_z_trunc = acf_z_avg[:idx_max]

    # Integración de la ACF
    Int_x = simps(acf_x_trunc, time_trunc)
    Int_y = simps(acf_y_trunc, time_trunc)
    Int_z = simps(acf_z_trunc, time_trunc)
    Int_avg = (Int_x + Int_y + Int_z) / 3.0

    # Green-Kubo: kappa = (V / (3 k_B T^2)) * Int
    GK_prefactor = V / (3.0 * kB_eV_K * T**2)
    kappa_x_SI = GK_prefactor * Int_x * FACTOR_CONVERSION
    kappa_y_SI = GK_prefactor * Int_y * FACTOR_CONVERSION
    kappa_z_SI = GK_prefactor * Int_z * FACTOR_CONVERSION
    kappa_avg_SI = GK_prefactor * Int_avg * FACTOR_CONVERSION

    # Resultados
    print("\n--- Resultados del Cálculo de Green-Kubo ---")
    print(f"Temperatura (T): {T:.2f} K")
    print(f"Volumen de Simulación (V): {V:.2f} Å^3")
    print(f"Tiempo de Truncamiento (tau_max): {tau_max:.1f} ps")
    print("-" * 40)
    print(f"Integral de ACF (X): {Int_x:.2e} [eV^2 Å^4 / ps]")
    print(f"Integral de ACF (Y): {Int_y:.2e} [eV^2 Å^4 / ps]")
    print(f"Integral de ACF (Z): {Int_z:.2e} [eV^2 Å^4 / ps]")
    print("-" * 40)
    print(f"Conductividad Térmica (W/(m K)):")
    print(f"kappa_x: {kappa_x_SI:.3f}")
    print(f"kappa_y: {kappa_y_SI:.3f}")
    print(f"kappa_z: {kappa_z_SI:.3f}")
    print(f"kappa_promedio: {kappa_avg_SI:.3f} 🌡️")

    # Graficar ACF
    plt.figure(figsize=(10, 6))
    plt.plot(time_trunc, acf_x_trunc, label=r'$\langle J_x(t) J_x(0) \rangle$')
    plt.plot(time_trunc, acf_y_trunc, label=r'$\langle J_y(t) J_y(0) \rangle$')
    plt.plot(time_trunc, acf_z_trunc, label=r'$\langle J_z(t) J_z(0) \rangle$')
    plt.axhline(0, color='grey', linestyle='--')
    plt.title('Función de Autocorrelación del Flujo de Calor')
    plt.xlabel('Tiempo (ps)')
    plt.ylabel('ACF [eV^2 Å^4 / ps^2]')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Graficar integral acumulativa
    cum_int_x = np.cumsum(acf_x_trunc) * dt_lammps
    cum_int_y = np.cumsum(acf_y_trunc) * dt_lammps
    cum_int_z = np.cumsum(acf_z_trunc) * dt_lammps
    plt.figure(figsize=(10, 6))
    plt.plot(time_trunc, cum_int_x, label='Integral Acumulativa X')
    plt.plot(time_trunc, cum_int_y, label='Integral Acumulativa Y')
    plt.plot(time_trunc, cum_int_z, label='Integral Acumulativa Z')
    plt.title('Integral Acumulativa de ACF para Convergencia')
    plt.xlabel('Tiempo (ps)')
    plt.ylabel('Integral [eV^2 Å^4 / ps]')
    plt.legend()
    plt.grid(True)
    plt.show()

    return {
        'kappa_x': kappa_x_SI,
        'kappa_y': kappa_y_SI,
        'kappa_z': kappa_z_SI,
        'kappa_avg': kappa_avg_SI,
        'Int_x': Int_x,
        'Int_y': Int_y,
        'Int_z': Int_z,
        'time': time_trunc,
        'acf_x': acf_x_trunc,
        'acf_y': acf_y_trunc,
        'acf_z': acf_z_trunc
    }

if __name__ == "__main__":
    try:
        sim_params = read_lammps_params()
        calculate_kappa(sim_params)
        sim_params['tau_max'] = 100.0  # Truncation time (ps)
        input_files = ['J0Jt_12345.dat', 'J0Jt_23456.dat', 'J0Jt_34567.dat']
        results = calculate_kappa_promedio(sim_params, input_files)
    except Exception as e:
        print(f"Fallo en el post-procesamiento: {e}")
