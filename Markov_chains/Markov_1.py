#!/usr/bin/env python
import random

def crear_matriz_transicion(estados, probabilidad_transicion):
    """
    Crea una matriz de transición para la cadena de Markov.
    estados: lista de símbolos posibles (bytes en este caso)
    probabilidad_transicion: dict con claves (estado_origen) y valores dict de prob. a estados destino
    """
    matriz = {}
    for estado in estados:
        # Obtenemos la dict con probabilidades para ese estado
        transiciones = probabilidad_transicion.get(estado, {})
        # Normalizamos las probabilidades para que sumen 1 (por si acaso)
        total = sum(transiciones.values())
        if total == 0:
            # Si no hay transiciones definidas, hacemos uniforme entre todos los estados
            matriz[estado] = {e: 1/len(estados) for e in estados}
        else:
            matriz[estado] = {e: transiciones.get(e, 0)/total for e in estados}
    return matriz

def siguiente_estado(estado_actual, matriz):
    """
    Dado el estado actual y la matriz de transición, devuelve el siguiente estado siguiendo probabilidades.
    """
    transiciones = matriz[estado_actual]
    estados = list(transiciones.keys())
    probabilidades = list(transiciones.values())
    return random.choices(estados, probabilidades)[0]

def generar_cadena_markov(matriz, estado_inicial, longitud):
    """
    Genera una secuencia de bytes siguiendo una cadena de Markov dada la matriz.
    """
    secuencia = [estado_inicial]
    estado_actual = estado_inicial
    for _ in range(longitud - 1):
        estado_actual = siguiente_estado(estado_actual, matriz)
        secuencia.append(estado_actual)
    return bytes(secuencia)

def guardar_archivo(nombre_archivo, datos):
    """
    Guarda los datos generados en un archivo binario.
    """
    with open(nombre_archivo, 'wb') as archivo:
        archivo.write(datos)

def main():
    # Definición del alfabeto de bytes (por simplicidad, bytes del 0x00 al 0x0F = 16 símbolos)
    estados = list(range(16))

    # Definición de una matriz de transiciones con probabilidad baja de entropía (más probabilidad en ciertos patrones)
    probabilidad_transicion = {
        0: {0: 0.7, 1: 0.1, 2: 0.1, 3: 0.1},
        1: {0: 0.6, 1: 0.2, 2: 0.1, 3: 0.1},
        2: {2: 0.8, 3: 0.1, 4: 0.1},
        3: {3: 0.7, 4: 0.3},
        # Estados restantes con transiciones Uniformes para evitar bloqueos
    }

    # Rellenar transiciones para estados no definidos con uniformes
    for estado in estados:
        if estado not in probabilidad_transicion:
            probabilidad_transicion[estado] = {e: 1/len(estados) for e in estados}

    matriz = crear_matriz_transicion(estados, probabilidad_transicion)

    # Estado inicial aleatorio
    estado_inicial = random.choice(estados)

    # Generar secuencia de 1 millón de bytes (puedes ajustar este tamaño)
    longitud = 1_000_000

    datos_generados = generar_cadena_markov(matriz, estado_inicial, longitud)

    # Guardar archivo resultado
    nombre_archivo = "archivo_ruido_markov.bin"
    guardar_archivo(nombre_archivo, datos_generados)

    print(f"Archivo generado: {nombre_archivo}")

if __name__ == "__main__":
    main()
