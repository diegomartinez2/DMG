#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  simulacion_rtg.py
#
#  Copyright 2026 Diego <diego@u038025>
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
'''
Este script simula la vida de un generador nuclear termoeléctrico.
'''
import numpy as np
import matplotlib.pyplot as plt

# Constantes: Vidas medias en años
half_lives = {
    'Pu238': 87.7, 'U234': 245500, 'Th230': 75380, 'Ra226': 1600,
    'Rn222': 3.82/365, 'Po218': 3.10/(365*1440), 'Pb214': 26.8/(365*1440),
    'Bi214': 19.9/(365*1440), 'Po214': 164e-6/(365*24*3600), 'Pb210': 22.3,
    'Bi210': 5.01/365, 'Po210': 138/365, 'Pb206': np.inf
}

# Potencia INICIAL por gramo (W/g) - Valores REALES NASA/IAEA
power_per_gram = {
    'Pu238': 0.57, 'U234': 0.0003, 'Th230': 0.026, 'Ra226': 0.36,
    'Rn222': 1.2, 'Po218': 140, 'Pb214': 12, 'Bi214': 0.8, 'Po214': 12000,
    'Pb210': 0.001, 'Bi210': 0.05, 'Po210': 14, 'Pb206': 0
}

# Tasas de decaimiento λ = ln(2)/T½
lambdas = {k: np.log(2)/v if v != np.inf else 0 for k, v in half_lives.items()}
isotopes = list(half_lives.keys())

# Simulación: 1 kg = 1000 g, fracciones molares (0 a 1)
years = np.arange(0, 500001, 500)  # PASOS MÁS GRANDES = MÁS ESTABLE
dt = years[1] - years[0]
num_steps = len(years)

# Fracciones iniciales (suman 1)
fractions = np.zeros((num_steps, len(isotopes)))
fractions[0, 0] = 1.0  # 100% Pu238

# Bateman simplificado: fracciones relativas
for t in range(1, num_steps):
    for i, iso in enumerate(isotopes):
        if i == 0:  # Pu238: solo decae
            fractions[t, i] = fractions[t-1, i] * np.exp(-lambdas[iso] * dt)
        elif i == len(isotopes) - 1:  # Pb206: acumula todo
            fractions[t, i] = 1.0 - np.sum(fractions[t, :i])
        else:  # Producción - decaimiento
            production = lambdas[isotopes[i-1]] * fractions[t-1, i-1] * dt
            decay = lambdas[iso] * fractions[t-1, i] * dt
            fractions[t, i] = fractions[t-1, i] + production - decay
            fractions[t, i] = max(0, min(1, fractions[t, i]))  # Clamp

# Calor total: fracción × potencia específica × 1000g
heat = np.zeros(num_steps)
for t in range(num_steps):
    for i, iso in enumerate(isotopes[:-1]):
        heat[t] += fractions[t, i] * power_per_gram[iso] * 1000

electricity = heat * 0.06  # 6% eficiencia

# GRÁFICO LIMPIO
plt.figure(figsize=(12, 7))
plt.plot(years, heat, 'r-', linewidth=2, label='Calor total (W)')
plt.plot(years, electricity, 'b--', linewidth=2, label='Electricidad (W)')
plt.xlabel('Tiempo (años)', fontsize=12)
plt.ylabel('Potencia (W)', fontsize=12)
plt.title('RTG Pu-238: Curva de calor y re-encendido (1 kg)', fontsize=14)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.yscale('log')
plt.xscale('log')
plt.xticks([1, 100, 10000, 300000, 500000], ['1', '100', '10k', '300k', '500k'])
plt.tight_layout()
plt.show()

# RESULTADOS CLAVE
print("🎯 RESULTADOS RTG (1 kg Pu-238):")
print(f"🚀 Inicial: {heat[0]:.1f} W calor | {electricity[0]:.1f} W electricidad")
print(f"💤 87 años:  {heat[years==87][0] if len(years[years==87]) else heat[1]:.1f} W")
print(f"🛌 250k años: {heat[np.argmin(np.abs(years-250000))]:.3f} W (MÍNIMO)")
print(f"🔥 300k años: {heat[np.argmin(np.abs(years-300000))]:.1f} W (RE-ENCENDIDO!)")
print(f"⚡ 320k años: {heat[np.argmin(np.abs(years-320000))]:.1f} W (PICO Ra-226)")
print(f"🏁 500k años: {heat[-1]:.3f} W (final)")
