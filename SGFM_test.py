#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
def surface_gf(omega, epsilon, t, eta=1e-5, max_iter=100, tol=1e-10):
    g = 1.0 / (omega + 1j * eta - epsilon) # Initial guess
    for _ in range(max_iter):
            g_new = 1.0 / (omega + 1j * eta - epsilon - t**2 * g)
            if np.abs(g_new - g) < tol:
                break
            g = g_new
    return g
# Parámetros
epsilon_A = 0.0
t_A = 1.0
epsilon_B = 0.5
t_B = 0.8
t_coupling = 0.2
omega = np.linspace(-3, 3, 500)
gf_interfaz = []
for w in omega:
    # Green superficial para cada material
    gA = surface_gf(w, epsilon_A, t_A)
    gB = surface_gf(w, epsilon_B, t_B)
    # Autoenergías
    Sigma_A = t_coupling**2 * gA
    Sigma_B = t_coupling**2 * gB
    # Green en la interfaz
    G00 = 1.0 / (w + 1j*1e-5 - 0.0 - Sigma_A - Sigma_B)
    gf_interfaz.append(-np.imag(G00) / np.pi)
plt.plot(omega, gf_interfaz)
plt.xlabel('$\omega$')
plt.ylabel('Densidad de estados (interfaz)')
plt.title('Densidad de estados en la interfaz (SGFM)')
plt.show()
