#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Algotitmo_Genético_para_XOR.py
#
#  Copyright 2025 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
# ---------------------------
# Importación de los módulos
# ---------------------------import numpy as np
import random
from typing import List, Tuple, Dict

# --- 1. Definición del Problema (Compuerta XOR) ---
# XOR es ideal porque no es linealmente separable.

# (Entrada, Salida Esperada)
XOR_DATA = [
    ([0, 0], [0]),
    ([0, 1], [1]),
    ([1, 0], [1]),
    ([1, 1], [0]),
]

# --- 2. La Red Neuronal (El Individuo) ---

class NeuralNetwork:
    """
    Una red neuronal pequeña con una capa oculta.
    Sus pesos y sesgos son el 'genoma' que será optimizado.
    """
    def __init__(self, input_size: int, hidden_size: int, output_size: int):
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.output_size = output_size

        # El número total de parámetros (genes) es la longitud del genoma
        # W1 (input -> hidden): input_size * hidden_size
        # B1 (bias hidden): hidden_size
        # W2 (hidden -> output): hidden_size * output_size
        # B2 (bias output): output_size
        self.genome_length = (input_size * hidden_size +
                              hidden_size +
                              hidden_size * output_size +
                              output_size)

        # Inicializar pesos aleatorios
        self.weights = np.random.randn(self.genome_length)

    def set_weights(self, new_weights: np.ndarray):
        """Asigna un nuevo genoma (vector de pesos) al individuo."""
        if new_weights.size != self.genome_length:
            raise ValueError(f"El genoma debe tener longitud {self.genome_length}, no {new_weights.size}")
        self.weights = new_weights

    def forward(self, input_data: np.ndarray) -> np.ndarray:
        """Propagación hacia adelante con función de activación sigmoide."""

        # Desempaquetar el vector de pesos (el genoma)
        i, h, o = self.input_size, self.hidden_size, self.output_size

        # W1: Pesos de entrada a capa oculta
        W1_end = i * h
        W1 = self.weights[0:W1_end].reshape(i, h)

        # B1: Sesgos de capa oculta
        B1_end = W1_end + h
        B1 = self.weights[W1_end:B1_end].reshape(1, h)

        # W2: Pesos de capa oculta a salida
        W2_end = B1_end + h * o
        W2 = self.weights[B1_end:W2_end].reshape(h, o)

        # B2: Sesgos de salida
        B2 = self.weights[W2_end:].reshape(1, o)

        # Capa oculta: H = sigmoid(Input @ W1 + B1)
        hidden_input = input_data @ W1 + B1
        hidden_output = 1 / (1 + np.exp(-hidden_input)) # Sigmoid

        # Capa de salida: O = sigmoid(H @ W2 + B2)
        output_input = hidden_output @ W2 + B2
        final_output = 1 / (1 + np.exp(-output_input)) # Sigmoid

        return final_output

# --- 3. El Algoritmo Genético (El Entorno Evolutivo) ---

def calculate_fitness(network: NeuralNetwork, data: List[Tuple[List[int], List[int]]]) -> float:
    """
    Evalúa la aptitud de la red. Una aptitud mayor es mejor.
    Usamos el error cuadrático medio negativo (negativo de la pérdida).
    """
    total_error = 0.0
    for input_val, expected_output in data:
        input_np = np.array(input_val).reshape(1, -1)
        expected_np = np.array(expected_output).reshape(1, -1)

        actual_output = network.forward(input_np)

        # Pérdida: Error cuadrático medio (MSE)
        error = np.mean((expected_np - actual_output)**2)
        total_error += error

    # El fitness es el NEGATIVO del error promedio (porque queremos MAXIMIZAR el fitness)
    # Un fitness de 0.0 es perfecto.
    # Usaremos 1.0 / (1.0 + total_error) para que el fitness esté entre 0 y 1.
    return 1.0 / (1.0 + total_error)

def select_parents(population: List[NeuralNetwork], fitnesses: List[float], num_parents: int) -> List[np.ndarray]:
    """
    Selección de padres por ranking: seleccionamos a los mejores individuos.
    """
    # Combinar fitness y genomas
    combined = sorted(zip(fitnesses, [n.weights for n in population]),
                      key=lambda x: x[0], reverse=True)

    # Retornar los genomas de los mejores 'num_parents'
    return [genome for fitness, genome in combined[:num_parents]]

def crossover(parent1: np.ndarray, parent2: np.ndarray) -> np.ndarray:
    """
    Cruce simple (de un punto). Combina genes de dos padres.
    """
    crossover_point = random.randint(1, len(parent1) - 1)

    child = np.concatenate([
        parent1[:crossover_point],
        parent2[crossover_point:]
    ])
    return child

def mutate(genome: np.ndarray, mutation_rate: float, mutation_strength: float) -> np.ndarray:
    """
    Aplica una pequeña perturbación aleatoria a los pesos con una probabilidad dada.
    """
    mutated_genome = np.copy(genome)

    for i in range(len(mutated_genome)):
        if random.random() < mutation_rate:
            # Añadir ruido Gaussiano (normal)
            mutated_genome[i] += random.gauss(0, mutation_strength)

    return mutated_genome

def run_neuroevolution(generations: int, pop_size: int,
                       mutation_rate: float, mutation_strength: float):
    """
    Bucle principal de la simulación evolutiva.
    """
    # Parámetros de la Red
    INPUT_SIZE = 2
    HIDDEN_SIZE = 4 # Pequeña capa oculta
    OUTPUT_SIZE = 1

    # Inicializar la población
    print(f"Genoma de longitud: {NeuralNetwork(INPUT_SIZE, HIDDEN_SIZE, OUTPUT_SIZE).genome_length}")
    population: List[NeuralNetwork] = [
        NeuralNetwork(INPUT_SIZE, HIDDEN_SIZE, OUTPUT_SIZE) for _ in range(pop_size)
    ]

    best_network = None
    best_fitness_history = []

    for gen in range(1, generations + 1):
        # 1. Evaluación de Aptitud (Fitness)
        fitnesses = [calculate_fitness(net, XOR_DATA) for net in population]

        max_fitness = max(fitnesses)
        best_index = fitnesses.index(max_fitness)
        best_network = population[best_index]
        best_fitness_history.append(max_fitness)

        if gen % 10 == 0 or gen == generations:
            print(f"Generación {gen:03d}/{generations} | Mejor Fitness: {max_fitness:.6f}")
            if max_fitness > 0.99: # Criterio de parada
                 print("¡Solución encontrada!")
                 break

        # 2. Selección de Padres
        # Seleccionar la mitad superior como padres
        num_parents = pop_size // 2
        parent_genomes = select_parents(population, fitnesses, num_parents)

        new_population_genomes = []

        # 3. Generación de Descendencia (Cruce y Mutación)
        # El bucle genera (pop_size - num_parents) nuevos hijos
        for _ in range(pop_size - num_parents):
            # Seleccionar dos padres al azar del pool de mejores padres
            p1_genome = random.choice(parent_genomes)
            p2_genome = random.choice(parent_genomes)

            # Cruce
            child_genome = crossover(p1_genome, p2_genome)

            # Mutación
            mutated_child_genome = mutate(child_genome, mutation_rate, mutation_strength)
            new_population_genomes.append(mutated_child_genome)

        # 4. Creación de la Nueva Población (Elitismo simple)
        # La nueva población incluye a los mejores padres (elitismo) y la nueva descendencia
        new_population = []
        for genome in parent_genomes:
            net = NeuralNetwork(INPUT_SIZE, HIDDEN_SIZE, OUTPUT_SIZE)
            net.set_weights(genome)
            new_population.append(net)

        for genome in new_population_genomes:
            net = NeuralNetwork(INPUT_SIZE, HIDDEN_SIZE, OUTPUT_SIZE)
            net.set_weights(genome)
            new_population.append(net)

        population = new_population

    # --- 5. Resultados Finales ---
    print("\n--- Resultados de la Mejor Red ---")

    if best_network:
        for input_val, expected in XOR_DATA:
            input_np = np.array(input_val).reshape(1, -1)
            output = best_network.forward(input_np)[0]

            # Redondear la salida (0 o 1) para comparación
            predicted = np.round(output)
            print(f"Entrada: {input_val} | Salida Esperada: {expected[0]} | Salida NN: {output[0]:.4f} (Predicción: {int(predicted[0])})")

    print(f"\nEntrenamiento completado en {gen} generaciones.")
    print("El fitness se calcula como 1 / (1 + Error Cuadrático Medio Total), donde 1.0 es perfecto.")


if __name__ == "__main__":
    # --- Parámetros del Algoritmo Genético ---
    NUM_GENERATIONS = 200     # Número máximo de generaciones a ejecutar
    POPULATION_SIZE = 50      # Tamaño de la población de redes
    MUTATION_RATE = 0.05      # Probabilidad de mutación por peso (similar a tu TASA_BASE_MUTACION)
    MUTATION_STRENGTH = 0.5   # Fuerza de la perturbación Gaussiana (similar a tu MAX_DESVIACION_MUTACION)

    # Ejecutar la evolución
    run_neuroevolution(
        generations=NUM_GENERATIONS,
        pop_size=POPULATION_SIZE,
        mutation_rate=MUTATION_RATE,
        mutation_strength=MUTATION_STRENGTH
    )
