#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
'''
For learning reasons, for real use is better to use the ones implemented in math libraries.
'''
# ==========================================
# 1. ALGORITMOS DE ORDENACIÓN (GENERADORES)
# ==========================================
# Se implementan como generadores (yield) para
# poder capturar cada estado intermedio del array.

def bubble_sort(arr):
    n = len(arr)
    for i in range(n):
        for j in range(0, n - i - 1):
            yield arr, j, j + 1  # Retorna el estado actual y los índices comparados
            if arr[j] > arr[j + 1]:
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
                yield arr, j, j + 1

def selection_sort(arr):
    n = len(arr)
    for i in range(n):
        min_idx = i
        for j in range(i + 1, n):
            yield arr, j, min_idx
            if arr[j] < arr[min_idx]:
                min_idx = j
        arr[i], arr[min_idx] = arr[min_idx], arr[i]
        yield arr, i, min_idx

def insertion_sort(arr):
    for i in range(1, len(arr)):
        key = arr[i]
        j = i - 1
        while j >= 0 and arr[j] > key:
            yield arr, j, j + 1
            arr[j + 1] = arr[j]
            j -= 1
            yield arr, j + 1, i
        arr[j + 1] = key
        yield arr, j + 1, i

def merge_sort(arr, start=0, end=None):
    if end is None:
        end = len(arr) - 1
    if start < end:
        mid = (start + end) // 2
        yield from merge_sort(arr, start, mid)
        yield from merge_sort(arr, mid + 1, end)

        # Proceso de fusión
        left = arr[start:mid + 1].copy()
        right = arr[mid + 1:end + 1].copy()
        i = j = 0
        for k in range(start, end + 1):
            yield arr, k, mid
            if i < len(left) and (j >= len(right) or left[i] <= right[j]):
                arr[k] = left[i]
                i += 1
            else:
                arr[k] = right[j]
                j += 1
            yield arr, k, mid

def quick_sort(arr, low=0, high=None):
    '''
    Note: non 'stable'
    '''
    if high is None:
        high = len(arr) - 1

    if low < high:
        # Partición de Lomuto
        pivot = arr[high]
        i = low - 1
        for j in range(low, high):
            yield arr, j, high
            if arr[j] < pivot:
                i += 1
                arr[i], arr[j] = arr[j], arr[i]
                yield arr, i, j
        arr[i + 1], arr[high] = arr[high], arr[i + 1]
        yield arr, i + 1, high
        p = i + 1

        yield from quick_sort(arr, low, p - 1)
        yield from quick_sort(arr, p + 1, high)

# ==========================================
# 2. COMPARACIÓN DE TIEMPOS DE EJECUCIÓN
# ==========================================

def medir_tiempos(algoritmos, tamaño_arr=200):
    """Mide el tiempo real que tarda cada algoritmo en ordenar un array idéntico."""
    datos_base = [random.randint(1, 1000) for _ in range(tamaño_arr)]
    tiempos = {}

    for nombre, (func, _) in algoritmos.items():
        arr_copia = datos_base.copy()
        inicio = time.perf_counter()
        # Consumir el generador completamente para ejecutar el ordenamiento
        for _ in func(arr_copia):
            pass
        fin = time.perf_counter()
        tiempos[nombre] = (fin - inicio) * 1000  # Convertir a milisegundos

    return tiempos

# ==========================================
# 3. VISUALIZACIÓN Y ANIMACIÓN
# ==========================================

def visualizar_algoritmo(nombre, func, datos):
    """Muestra la animación paso a paso de un algoritmo individual."""
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.canvas.manager.set_window_title(f"Visualizador: {nombre}")

    arr_copia = datos.copy()
    barras = ax.bar(range(len(arr_copia)), arr_copia, color="steelblue", align="edge")

    ax.set_title(f"Algoritmo: {nombre}", fontsize=14, fontweight="bold")
    ax.set_xlim(0, len(arr_copia))
    ax.set_ylim(0, max(arr_copia) * 1.1)

    texto_pasos = ax.text(0.02, 0.95, "", transform=ax.transAxes, fontsize=11)
    generador = func(arr_copia)
    pasos = [0]

    def update(frame):
        try:
            arr, idx1, idx2 = next(generador)
            pasos[0] += 1
            for idx, bar in enumerate(barras):
                bar.set_height(arr[idx])
                if idx in (idx1, idx2):
                    bar.set_color("crimson")  # Resaltar elementos siendo comparados
                else:
                    bar.set_color("steelblue")
            texto_pasos.set_text(f"Operaciones/Comparaciones: {pasos[0]}")
        except StopIteration:
            # Al finalizar, colorear todo de verde
            for bar in barras:
                bar.set_color("forestgreen")
            anim.event_source.stop()
        return barras

    anim = FuncAnimation(fig, update, frames=None, interval=15, repeat=False, cache_frame_data=False)
    plt.show()

def mostrar_comparativa_tiempos(tiempos):
    """Muestra un gráfico de barras comparando los tiempos de ejecución."""
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.canvas.manager.set_window_title("Comparativa de Tiempos")

    nombres = list(tiempos.keys())
    valores = list(tiempos.values())

    barras = ax.bar(nombres, valores, color=["#e74c3c", "#e67e22", "#f1c40f", "#2ecc71", "#3498db"])

    ax.set_ylabel("Tiempo de ejecución (ms)", fontsize=12)
    ax.set_title("Comparativa de Rendimiento (Array de 200 elementos)", fontsize=14, fontweight="bold")
    ax.grid(axis="y", linestyle="--", alpha=0.7)

    for bar in barras:
        altura = bar.get_height()
        ax.annotate(f'{altura:.2f} ms',
                    xy=(bar.get_x() + bar.get_width() / 2, altura),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontweight='bold')

    plt.tight_layout()
    plt.show()

# ==========================================
# 4. MENÚ PRINCIPAL
# ==========================================

def main():
    algoritmos = {
        "Bubble Sort": (bubble_sort, "$O(n^2)$"),
        "Selection Sort": (selection_sort, "$O(n^2)$"),
        "Insertion Sort": (insertion_sort, "$O(n^2)$"),
        "Merge Sort": (merge_sort, "$O(n \\log n)$"),
        "Quick Sort": (quick_sort, "$O(n \\log n)$")
    }

    TAMAÑO_ARRAY_ANIMACION = 40  # Tamaño reducido para apreciar la animación paso a paso
    datos_animacion = [random.randint(10, 100) for _ in range(TAMAÑO_ARRAY_ANIMACION)]

    while True:
        print("\n" + "="*45)
        print("    VISUALIZADOR DE ALGORITMOS DE ORDENACIÓN")
        print("="*45)
        for i, (nombre, (_, complejidad)) in enumerate(algoritmos.items(), 1):
            print(f"{i}. Animación: {nombre:<15} [Complejidad: {complejidad}]")
        print("6. Ver comparativa de tiempos de ejecución (Benchmarks)")
        print("7. Salir")

        opcion = input("\nSelecciona una opción (1-7): ").strip()

        if opcion in [str(i) for i in range(1, 6)]:
            idx = int(opcion) - 1
            nombre = list(algoritmos.keys())[idx]
            func, _ = algoritmos[nombre]
            print(f"\nAbriendo ventana de animación para {nombre}...")
            visualizar_algoritmo(nombre, func, datos_animacion)
        elif opcion == "6":
            print("\nCalculando tiempos de ejecución en un array de 200 elementos aleatorios...")
            tiempos = medir_tiempos(algoritmos, tamaño_arr=200)
            mostrar_comparativa_tiempos(tiempos)
        elif opcion == "7":
            print("\n¡Hasta luego!")
            break
        else:
            print("Opción inválida. Inténtalo de nuevo.")

if __name__ == "__main__":
    main()
