import numpy as np
import matplotlib.pyplot as plt

def calculate_efficiency(ZT, Th_celsius, Tc_celsius):
    # Convertir a Kelvin
    Th = Th_celsius + 273.15
    Tc = Tc_celsius + 273.15
    T_avg = (Th + Tc) / 2

    # Eficiencia de Carnot
    carnot_eff = (Th - Tc) / Th

    # Factor de mérito del material
    m = np.sqrt(1 + ZT)
    reduction_factor = (m - 1) / (m + (Tc / Th))

    return carnot_eff * reduction_factor

# Rango de ZT de 0 a 20 (el límite práctico vs teórico)
zt_range = np.linspace(0.01, 20, 500)

# Configuraciones de temperatura (Th, Tc) en Celsius
scenarios = [
    (500, 20, 'Industrial/Motor (500°C)'),
    (200, 20, 'Calor residual bajo (200°C)'),
    (60, 20, 'Electrónica/CPU (60°C)')
]

plt.figure(figsize=(10, 6))

for Th, Tc, label in scenarios:
    effs = [calculate_efficiency(zt, Th, Tc) * 100 for zt in zt_range]
    carnot = ((Th + 273.15) - (Tc + 273.15)) / (Th + 273.15) * 100

    line, = plt.plot(zt_range, effs, label=f'Eficiencia {label}', linewidth=2)
    plt.axhline(y=carnot, color=line.get_color(), linestyle='--', alpha=0.5,
                label=f'Límite Carnot ({label})')

# Formateo de la gráfica
plt.title('Eficiencia Termoeléctrica vs Factor de Mérito (ZT)', fontsize=14)
plt.xlabel('Cifra de mérito (ZT)', fontsize=12)
plt.ylabel('Eficiencia (%)', fontsize=12)
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlim(0, 20)
plt.ylim(0, 70)
plt.tight_layout()

print("Generando gráfica... Si tienes un entorno con soporte visual (como Jupyter o IDE local), verás la figura.")
plt.show()
