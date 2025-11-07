# -------------------------------------------------
# file: collatz_visualizer.py
# -------------------------------------------------
"""
Simulación visual de la conjetura de Collatz.
Muestra cómo diferentes números convergen al ciclo 4-2-1.
"""

import matplotlib.pyplot as plt
from typing import List, Tuple


def collatz_sequence(n: int) -> List[int]:
    """
    Genera la secuencia de Collatz completa hasta llegar a 1.
    Incluye el 1 final para cerrar el ciclo.
    """
    seq = [n]
    while n != 1:
        if n % 2 == 0:
            n //= 2
        else:
            n = 3 * n + 1
        seq.append(n)
    return seq


def plot_collatz_sequences(start_numbers: List[int], max_steps: int = 200) -> None:
    """
    Grafica múltiples secuencias de Collatz.
    - Cada línea es un número inicial.
    - Eje Y: valor del número (escala logarítmica).
    - Eje X: paso en la secuencia.
    """
    plt.figure(figsize=(12, 7))
    colors = plt.cm.tab10.colors

    for i, start in enumerate(start_numbers):
        seq = collatz_sequence(start)
        # Trunca si es muy larga
        if len(seq) > max_steps:
            seq = seq[:max_steps]
        steps = list(range(len(seq)))
        plt.plot(steps, seq, marker='o', markersize=3, linewidth=1.5,
                 color=colors[i % len(colors)], label=f'n={start}')

    # Resaltar el ciclo 4-2-1
    cycle = [4, 2, 1]
    plt.axhline(4, color='red', linestyle='--', alpha=0.6, label='Ciclo 4-2-1')
    plt.axhline(2, color='red', linestyle='--', alpha=0.6)
    plt.axhline(1, color='red', linestyle='--', alpha=0.6)

    plt.yscale('log')
    plt.xlabel('Paso en la secuencia')
    plt.ylabel('Valor (escala log)')
    plt.title('Conjetura de Collatz: Trayectorias hacia el ciclo 4→2→1')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


# -------------------------------------------------
# Ejemplo de uso
# -------------------------------------------------
if __name__ == "__main__":
    # Números interesantes: potencia de 2, número grande, número con pico alto
    starters = [6, 27, 100, 63728127]  # 27 tiene el pico más alto conocido
    plot_collatz_sequences(starters)
