import numpy as np

class AntColonyTSP:
    def __init__(self, coords, n_ants=15, n_iterations=100, alpha=1.0, beta=3.0, evaporation_rate=0.5, Q=100.0):
        """
        Inicializa el optimizador de colonia de hormigas para el TSP.

        :param coords: Lista de tuplas con las coordenadas (x, y) de las ciudades.
        :param n_ants: Número de hormigas por iteración.
        :param n_iterations: Número total de generaciones/iteraciones.
        :param alpha: Importancia relativa de la feromona (memoria colectiva).
        :param beta: Importancia relativa de la heurística (visibilidad/distancia inversa).
        :param evaporation_rate: Tasa de evaporación de la feromona (rho).
        :param Q: Constante de depósito de feromona.
        """
        self.coords = coords
        self.n_cities = len(coords)
        self.n_ants = n_ants
        self.n_iterations = n_iterations
        self.alpha = alpha
        self.beta = beta
        self.rho = evaporation_rate
        self.Q = Q

        # Matriz de distancias euclidianas entre todas las ciudades
        self.dist_matrix = np.zeros((self.n_cities, self.n_cities))
        for i in range(self.n_cities):
            for j in range(self.n_cities):
                self.dist_matrix[i][j] = np.linalg.norm(np.array(coords[i]) - np.array(coords[j]))

        # Matriz heurística (eta = 1 / distancia). Evitamos división por cero en la diagonal.
        with np.errstate(divide='ignore'):
            self.eta = 1.0 / self.dist_matrix
            np.fill_diagonal(self.eta, 0)

        # Inicialización de la matriz de feromonas con un valor pequeño constante
        self.pheromone = np.ones((self.n_cities, self.n_cities)) * 0.1

    def run(self):
        best_tour = None
        best_distance = float('inf')

        for iteration in range(self.n_iterations):
            all_tours = []
            all_distances = []

            # Cada hormiga construye un recorrido completo
            for _ in range(self.n_ants):
                tour = self._construct_tour()
                distance = self._calculate_tour_distance(tour)
                all_tours.append(tour)
                all_distances.append(distance)

                # Actualizar la mejor solución global encontrada
                if distance < best_distance:
                    best_distance = distance
                    best_tour = tour

            # Actualizar las feromonas en el grafo
            self._update_pheromones(all_tours, all_distances)

        return best_tour, best_distance

    def _construct_tour(self):
        start_node = np.random.randint(0, self.n_cities)
        tour = [start_node]
        unvisited = set(range(self.n_cities))
        unvisited.remove(start_node)

        current_node = start_node
        while unvisited:
            next_node = self._select_next_node(current_node, unvisited)
            tour.append(next_node)
            unvisited.remove(next_node)
            current_node = next_node

        return tour

    def _select_next_node(self, current_node, unvisited):
        unvisited_list = list(unvisited)

        # Obtener feromonas y heurísticas para los nodos no visitados
        pheromones = np.array([self.pheromone[current_node][j] for j in unvisited_list])
        heuristics = np.array([self.eta[current_node][j] for j in unvisited_list])

        # Calcular probabilidades según la regla de transición
        probabilities = (pheromones ** self.alpha) * (heuristics ** self.beta)
        sum_prob = np.sum(probabilities)

        if sum_prob == 0:
            return np.random.choice(unvisited_list)

        probabilities /= sum_prob

        # Seleccionar el siguiente nodo estocásticamente basado en las probabilidades
        return np.random.choice(unvisited_list, p=probabilities)

    def _calculate_tour_distance(self, tour):
        distance = 0
        for i in range(len(tour)):
            from_node = tour[i]
            to_node = tour[(i + 1) % len(tour)]  # Cierra el ciclo de vuelta al inicio
            distance += self.dist_matrix[from_node][to_node]
        return distance

    def _update_pheromones(self, all_tours, all_distances):
        # 1. Evaporación global de la feromona existente
        self.pheromone *= (1.0 - self.rho)

        # 2. Depósito de nueva feromona por cada hormiga (inversamente proporcional a la distancia)
        for tour, distance in zip(all_tours, all_distances):
            deposit = self.Q / distance
            for i in range(len(tour)):
                from_node = tour[i]
                to_node = tour[(i + 1) % len(tour)]
                self.pheromone[from_node][to_node] += deposit
                self.pheromone[to_node][from_node] += deposit  # Grafo simétrico

# --- Ejemplo de ejecución ---
if __name__ == "__main__":
    # Coordenadas aleatorias de 10 ciudades en un plano 2D
    np.random.seed(42)
    cities_coords = [(np.random.randint(0, 100), np.random.randint(0, 100)) for _ in range(10)]

    print("Coordenadas de las ciudades:")
    for idx, coord in enumerate(cities_coords):
        print(f"Ciudad {idx}: {coord}")
    print("-" * 40)

    # Configurar y ejecutar el algoritmo de la colonia de hormigas
    aco = AntColonyTSP(
        coords=cities_coords,
        n_ants=20,
        n_iterations=150,
        alpha=1.0,   # Peso de la feromona
        beta=4.0,    # Peso de la distancia (prioriza caminos cortos)
        evaporation_rate=0.5,
        Q=100.0
    )

    best_tour, best_distance = aco.run()

    print("¡Optimización completada!")
    print(f"Mejor ruta (orden de visita de ciudades): {best_tour}")
    print(f"Distancia total óptima: {round(best_distance, 2)}")
