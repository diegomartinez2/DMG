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

def main():
    """Función principal para ejecutar la simulación."""
    size = 10  # Tamaño de la cuadrícula
    steps = 50  # Número máximo de pasos
    grid = Grid(size, seed=42)  # Cuadrícula inicial
    evolution = Evolution(base_sigma=0.5, neighbor_weight=0.5)  # Evolución
    visualizer = Visualizer(graphical=False)  # Visualización en terminal (cambiar a True para gráfica)

    # Visualizar estado inicial
    visualizer.display(grid, 0)

    # Evolucionar y visualizar
    history = evolution.evolve(grid, steps)
    for step in range(1, len(history) + 1):
        visualizer.display(grid, step)

    # Mostrar evolución del valor medio absoluto
    plt.ioff()
    plt.figure()
    plt.plot(history)
    plt.xlabel('Paso')
    plt.ylabel('Valor medio absoluto')
    plt.title('Evolución de la entropía aproximada')
    plt.show()

if __name__ == "__main__":
    main()
