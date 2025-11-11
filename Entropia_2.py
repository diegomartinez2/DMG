import numpy as np
import matplotlib.pyplot as plt
import time

class Grid:
    """Gestiona la cuadrícula y sus valores."""
    def __init__(self, size, seed=None):
        if seed is not None:
            np.random.seed(seed)
        self.size = size
        self.values = np.random.uniform(-1, 1, (size, size))

    def get_neighbors(self, i, j):
        size = self.size
        neighbors = [
            self.values[(i-1) % size, j],
            self.values[(i+1) % size, j],
            self.values[i, (j-1) % size],
            self.values[i, (j+1) % size]
        ]
        return neighbors

class EntropyCalculator:
    """Calcula la entropía de Shannon de la cuadrícula."""
    def __init__(self, num_bins=20):
        self.num_bins = num_bins

    def calculate(self, grid):
        values = grid.values.flatten()
        hist, bin_edges = np.histogram(values, bins=self.num_bins, range=(-1, 1), density=True)
        probs = hist / np.sum(hist)
        probs = probs[probs > 0]
        entropy = -np.sum(probs * np.log2(probs)) if probs.size > 0 else 0
        return entropy

class Evolution:
    """Gestiona la evolución de la cuadrícula."""
    def __init__(self, base_sigma=0.5, neighbor_weight=0.5):
        self.base_sigma = base_sigma
        self.neighbor_weight = neighbor_weight
        self.entropy_calculator = EntropyCalculator(num_bins=20)

    def update_cell(self, grid, i, j):
        neighbors = grid.get_neighbors(i, j)
        neighbor_mean_abs = np.mean([abs(n) for n in neighbors])
        sigma = self.base_sigma * (1 - self.neighbor_weight * neighbor_mean_abs)
        sigma = max(sigma, 0.01)
        new_value = np.random.normal(0, sigma)
        return np.clip(new_value, -1, 1)

    def evolve(self, grid, steps=100, equilibrium_threshold=0.01):
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
            if step > 0 and abs(mean_history[-1] - mean_history[-2]) < equilibrium_threshold:
                print(f"Equilibrio alcanzado en el paso {step+1}")
                break
        return mean_history, entropy_history

class Visualizer:
    """Gestiona la visualización de la cuadrícula."""
    def __init__(self, graphical=False):
        self.graphical = graphical
        if self.graphical:
            plt.ion()

    def display(self, grid, step):
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
            time.sleep(0.1)

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
