#!/usr/bin/env python
#This script doesn't work, need fix
import numpy as np
import matplotlib.pyplot as plt
import time

class Grid:
    """Gestiona la cuadrícula y sus valores. Cumple SRP al manejar solo la estructura de datos."""
    def __init__(self, size, seed=None):
        """
        Inicializa una cuadrícula NxN con valores aleatorios entre -1 y 1.
        Args:
            size (int): Tamaño de la cuadrícula (NxN).
            seed (int): Semilla para reproducibilidad.
        """
        if seed is not None:
            np.random.seed(seed)
        self.size = size
        self.values = np.random.uniform(-1, 1, (size, size))

    def get_neighbors(self, i, j):
        """
        Devuelve los valores de las celdas adyacentes con condiciones periódicas.
        Args:
            i, j (int): Índices de la celda.
        Returns:
            list: Valores de las celdas adyacentes (arriba, abajo, izquierda, derecha).
        """
        size = self.size
        neighbors = [
            self.values[(i-1) % size, j],  # Arriba
            self.values[(i+1) % size, j],  # Abajo
            self.values[i, (j-1) % size],  # Izquierda
            self.values[i, (j+1) % size]   # Derecha
        ]
        return neighbors

class Evolution:
    """Gestiona la evolución de la cuadrícula. Cumple SRP al manejar solo las actualizaciones."""
    def __init__(self, base_sigma=0.5, neighbor_weight=0.5):
        """
        Inicializa parámetros de evolución.
        Args:
            base_sigma (float): Desviación estándar base para la gaussiana.
            neighbor_weight (float): Peso de los valores vecinos en el cálculo de sigma.
        """
        self.base_sigma = base_sigma
        self.neighbor_weight = neighbor_weight

    def update_cell(self, grid, i, j):
        """
        Actualiza una celda con una distribución gaussiana cuya sigma depende de los vecinos.
        Args:
            grid (Grid): Cuadrícula con los valores.
            i, j (int): Índices de la celda.
        Returns:
            float: Nuevo valor de la celda, limitado entre -1 y 1.
        """
        neighbors = grid.get_neighbors(i, j)
        # Calculamos sigma: más estrecha si los valores vecinos están cerca de 0
        neighbor_mean_abs = np.mean([abs(n) for n in neighbors])
        sigma = self.base_sigma * (1 - self.neighbor_weight * neighbor_mean_abs)
        sigma = max(sigma, 0.01)  # Evitar sigma demasiado pequeño
        new_value = np.random.normal(0, sigma)
        return np.clip(new_value, -1, 1)

    def evolve(self, grid, steps=100, equilibrium_threshold=0.01):
        """
        Evoluciona la cuadrícula durante varios pasos o hasta alcanzar el equilibrio.
        Args:
            grid (Grid): Cuadrícula a evolucionar.
            steps (int): Número máximo de pasos.
            equilibrium_threshold (float): Umbral para considerar equilibrio.
        Returns:
            list: Historia de los valores medios absolutos para análisis.
        """
        history = []
        for step in range(steps):
            old_values = grid.values.copy()
            for i in range(grid.size):
                for j in range(grid.size):
                    grid.values[i, j] = self.update_cell(grid, i, j)
            mean_abs_value = np.mean(np.abs(grid.values))
            history.append(mean_abs_value)
            # Criterio de equilibrio: cambio promedio pequeño
            if step > 0 and abs(history[-1] - history[-2]) < equilibrium_threshold:
                print(f"Equilibrio alcanzado en el paso {step+1}")
                break
        return history

class Visualizer:
    """Gestiona la visualización de la cuadrícula. Cumple SRP al manejar solo la representación."""
    def __init__(self, graphical=False):
        """
        Inicializa el visualizador.
        Args:
            graphical (bool): Si True, usa gráficos de colores; si False, usa terminal.
        """
        self.graphical = graphical
        if self.graphical:
            plt.ion()  # Modo interactivo para actualización en tiempo real

    def display(self, grid, step):
        """
        Muestra la cuadrícula en la terminal o gráficamente.
        Args:
            grid (Grid): Cuadrícula a visualizar.
            step (int): Paso actual para el título.
        """
        if self.graphical:
            plt.clf()
            plt.imshow(np.abs(grid.values), cmap='hot', vmin=0, vmax=1)
            plt.colorbar(label='Valor absoluto')
            plt.title(f'Paso {step}')
            plt.pause(0.1)
        else:
            print(f"\nPaso {step}:")
            for row in np.abs(grid.values):
                print(" ".join(f"{x:.2f}" for x in row))
            time.sleep(0.1)  # Pausa para visualización en terminal

class Evolution:
    def __init__(self, base_sigma=0.5, neighbor_weight=0.5):
        self.base_sigma = base_sigma
        self.neighbor_weight = neighbor_weight
        self.entropy_calculator = EntropyCalculator(num_bins=20)  # Añadimos el calculador

    def evolve(self, grid, steps=100, equilibrium_threshold=0.01):
        """
        Evoluciona la cuadrícula y calcula métricas (valor medio absoluto y entropía).
        Returns:
            tuple: (historia de valores medios absolutos, historia de entropía)
        """
        mean_history = []
        entropy_history = []
        for step in range(steps):
            old_values = grid.values.copy()
            for i in range(grid.size):
                for j in range(grid.size):
                    grid.values[i, j] = self.update_cell(grid, i, j)
            mean_abs_value = np.mean(np.abs(grid.values))
            entropy = self.entropy_calculator.calculate(grid)
            mean_history.append(mean_abs_value)
            entropy_history.append(entropy)
            # Criterio de equilibrio
            if step > 0 and abs(mean_history[-1] - mean_history[-2]) < equilibrium_threshold:
                print(f"Equilibrio alcanzado en el paso {step+1}")
                break
        return mean_history, entropy_history

class EntropyCalculator:
    """Calcula la entropía de Shannon de la cuadrícula. Cumple SRP al manejar solo el cálculo de entropía."""
    def __init__(self, num_bins=20):
        """
        Inicializa el calculador de entropía.
        Args:
            num_bins (int): Número de bins para discretizar los valores entre -1 y 1.
        """
        self.num_bins = num_bins

    def calculate(self, grid):
        """
        Calcula la entropía de Shannon de los valores de la cuadrícula.
        Args:
            grid (Grid): Cuadrícula con los valores.
        Returns:
            float: Entropía en bits (usando log base 2).
        """
        # Discretizar valores en bins
        values = grid.values.flatten()
        hist, bin_edges = np.histogram(values, bins=self.num_bins, range=(-1, 1), density=True)
        # Calcular probabilidades (normalizadas)
        probs = hist / np.sum(hist)
        # Evitar log(0) filtrando probabilidades nulas
        probs = probs[probs > 0]
        # Entropía de Shannon: -sum(p * log2(p))
        entropy = -np.sum(probs * np.log2(probs)) if probs.size > 0 else 0
        return entropy


def main():
    size = 10
    steps = 50
    grid = Grid(size, seed=42)
    evolution = Evolution(base_sigma=0.5, neighbor_weight=0.5)
    visualizer = Visualizer(graphical=False)

    visualizer.display(grid, 0)
    mean_history, entropy_history = evolution.evolve(grid, steps)
    for step in range(1, len(mean_history) + 1):
        visualizer.display(grid, step)

    # Graficar evolución
    plt.ioff()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
    ax1.plot(mean_history)
    ax1.set_xlabel('Paso')
    ax1.set_ylabel('Valor medio absoluto')
    ax1.set_title('Evolución del valor medio absoluto')
    ax2.plot(entropy_history)
    ax2.set_xlabel('Paso')
    ax2.set_ylabel('Entropía (bits)')
    ax2.set_title('Evolución de la entropía de Shannon')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
