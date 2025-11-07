# -------------------------------------------------
# file: optimal_transport_demo.py
# -------------------------------------------------
"""
Demostración de Transporte Óptimo:
1. Clásico (1D discreto)
2. Con restricción de entrelazamiento (simulación simple)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import linprog
from typing import Tuple


# -------------------------------------------------
# 1. Transporte Óptimo Clásico
# -------------------------------------------------
def optimal_transport_classical(
    mu: np.ndarray, nu: np.ndarray, cost_matrix: np.ndarray
) -> Tuple[np.ndarray, float]:
    """
    Resuelve transporte óptimo discreto con programación lineal.
    - mu: distribución inicial (suma a 1)
    - nu: distribución final (suma a 1)
    - cost_matrix: c[i,j] = costo de mover de i a j
    Retorna: (matriz de transporte, costo total)
    """
    n, m = cost_matrix.shape
    # Aplanar para LP
    c = cost_matrix.ravel()

    # Restricciones: suma filas = mu, suma columnas = nu
    A_eq = []
    b_eq = []
    # Filas
    for i in range(n):
        row = np.zeros((n, m))
        row[i, :] = 1
        A_eq.append(row.ravel())
        b_eq.append(mu[i])
    # Columnas
    for j in range(m):
        col = np.zeros((n, m))
        col[:, j] = 1
        A_eq.append(col.ravel())
        b_eq.append(nu[j])

    A_eq = np.array(A_eq)
    b_eq = np.array(b_eq)

    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method='highs')
    if not res.success:
        raise ValueError("No se pudo resolver el transporte óptimo")

    transport = res.x.reshape(n, m)
    total_cost = res.fun
    return transport, total_cost


# -------------------------------------------------
# 2. Transporte con Restricción de Entrelazamiento
# -------------------------------------------------
def optimal_transport_quantum(
    mu: np.ndarray, nu: np.ndarray, cost_matrix: np.ndarray, max_corr: float = 0.9
) -> Tuple[np.ndarray, float]:
    """
    Añade restricción: correlación entre pares no puede exceder max_corr.
    Simula entrelazamiento: si mueves mucho de i a j, no puedes mover de i' a j'.
    """
    n, m = cost_matrix.shape
    c = cost_matrix.ravel()

    A_eq, b_eq = [], []
    # Igual que clásico
    for i in range(n):
        row = np.zeros((n, m))
        row[i, :] = 1
        A_eq.append(row.ravel())
        b_eq.append(mu[i])
    for j in range(m):
        col = np.zeros((n, m))
        col[:, j] = 1
        A_eq.append(col.ravel())
        b_eq.append(nu[j])

    # Restricción de entrelazamiento: |P[i,j] - P[i',j']| <= max_corr para pares entrelazados
    A_ub, b_ub = [], []
    entangled_pairs = [(0, 1)]  # Ej: partículas 0 y 1 están entrelazadas
    for i1, i2 in entangled_pairs:
        for j1 in range(m):
            for j2 in range(m):
                if j1 == j2:
                    continue
                row = np.zeros(n * m)
                row[i1 * m + j1] = 1
                row[i2 * m + j2] = -1
                A_ub.append(row)
                b_ub.append(max_corr)
                A_ub.append(-row)
                b_ub.append(max_corr)

    res = linprog(
        c, A_eq=A_eq, b_eq=b_eq, A_ub=A_ub, b_ub=b_ub,
        bounds=(0, None), method='highs'
    )
    if not res.success:
        raise ValueError("Transporte cuántico falló")

    transport = res.x.reshape(n, m)
    return transport, res.fun


# -------------------------------------------------
# Visualización
# -------------------------------------------------
def plot_transport(transport: np.ndarray, title: str):
    plt.figure(figsize=(6, 5))
    im = plt.imshow(transport, cmap='Blues', origin='lower')
    plt.colorbar(im, label='Masa transportada')
    plt.title(title)
    plt.xlabel('Destino j')
    plt.ylabel('Origen i')
    plt.xticks(range(transport.shape[1]))
    plt.yticks(range(transport.shape[0]))
    for i in range(transport.shape[0]):
        for j in range(transport.shape[1]):
            plt.text(j, i, f'{transport[i,j]:.2f}', ha='center', va='center')
    plt.tight_layout()
    plt.show()


# -------------------------------------------------
# Ejemplo completo
# -------------------------------------------------
if __name__ == "__main__":
    # Distribuciones: 4 puntos origen, 3 destino
    mu = np.array([0.4, 0.3, 0.2, 0.1])  # suma = 1
    nu = np.array([0.0, 0.5, 0.5])      # suma = 1

    # Costo: distancia euclidiana en 1D
    x_orig = np.array([0, 1, 3, 6])
    x_dest = np.array([2, 4, 5])
    cost = np.abs(x_orig[:, None] - x_dest[None, :])

    print("Transporte Óptimo Clásico")
    T_classic, cost_classic = optimal_transport_classical(mu, nu, cost)
    print(f"Costo total: {cost_classic:.3f}")
    plot_transport(T_classic, "Transporte Óptimo Clásico")

    print("\nTransporte con Restricción Cuántica (entrelazamiento 0-1)")
    T_quantum, cost_quantum = optimal_transport_quantum(mu, nu, cost, max_corr=0.7)
    print(f"Costo total (con restricción): {cost_quantum:.3f}")
    plot_transport(T_quantum, "Transporte con Límite de Entrelazamiento")
