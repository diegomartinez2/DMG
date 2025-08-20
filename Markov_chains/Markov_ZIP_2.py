import random

def crear_matriz_transicion(estados, probabilidad_transicion):
    """
    Crea una matriz de transición para la cadena de Markov.
    estados: lista de símbolos posibles (bytes en este caso)
    probabilidad_transicion: dict con claves (estado_origen) y valores dict de prob. a estados destino
    """
    matriz = {}
    for estado in estados:
        transiciones = probabilidad_transicion.get(estado, {})
        total = sum(transiciones.values())
        if total == 0:
            # Transición uniforme si no definido
            matriz[estado] = {e: 1/len(estados) for e in estados}
        else:
            matriz[estado] = {e: transiciones.get(e, 0)/total for e in estados}
    return matriz

def siguiente_estado(estado_actual, matriz):
    """
    Devuelve el siguiente estado según la matriz de transición y estado actual.
    """
    transiciones = matriz[estado_actual]
    estados = list(transiciones.keys())
    probabilidades = list(transiciones.values())
    return random.choices(estados, probabilidades)[0]

def generar_cadena_markov(matriz, estado_inicial, longitud):
    """
    Genera secuencia de bytes siguiendo la cadena de Markov.
    """
    secuencia = [estado_inicial]
    estado_actual = estado_inicial
    for _ in range(longitud - 1):
        estado_actual = siguiente_estado(estado_actual, matriz)
        secuencia.append(estado_actual)
    return bytes(secuencia)

def guardar_archivo(nombre_archivo, datos):
    """
    Guarda bytes en archivo binario.
    """
    with open(nombre_archivo, 'wb') as archivo:
        archivo.write(datos)

def generar_encabezado_zip():
    """
    Encabezado típico ZIP (4 bytes de firma)
    """
    return bytes([0x50, 0x4B, 0x03, 0x04])

def main():
    # Rango completo bytes (0x00 a 0xFF)
    estados = list(range(256))

    # Ejemplo simplificado: probabilidades que favorecen la aparición de los bytes de "PK\x03\x04"
    # y algunos rangos comunes en binarios (por simplicidad se priorizan sólo algunos bytes)
    probabilidad_transicion = {
        0x50: {0x4B: 0.6, 0x50: 0.2, 0x00: 0.2},  # P(0x50) muy probable seguir con 0x4B (K)
        0x4B: {0x03: 0.7, 0x50: 0.1, 0x4B: 0.2},
        0x03: {0x04: 0.7, 0x00: 0.3},
        0x04: {0x50: 0.4, 0x00: 0.6},
    }
    # Para el resto: distribución uniforme para evitar estados muertos
    for estado in estados:
        if estado not in probabilidad_transicion:
            # Algunos bytes con mayor probabilidad para imitar rango común en binarios
            # Ejemplo: bytes de cero a 32 (control) y rango alto para bytes aleatorios
            transiciones = {}
            for e in estados:
                # Suavizar distribución, favorecer algunos bloques:
                if e in (0x00, 0x20, 0xFF, 0x7F):  # común en binarios o caracteres de control
                    transiciones[e] = 3.0
                else:
                    transiciones[e] = 1.0
            # Normalizar luego dentro de crear_matriz_transicion
            probabilidad_transicion[estado] = transiciones

    matriz = crear_matriz_transicion(estados, probabilidad_transicion)

    # Estado inicial que comienza la firma ZIP 'P' = 0x50
    estado_inicial = 0x50

    tamaño_objetivo = 154  # bytes totales del archivo falso
    encabezado = generar_encabezado_zip()
    longitud_ruido = tamaño_objetivo - len(encabezado)

    datos_ruido = generar_cadena_markov(matriz, estado_inicial, longitud_ruido)
    datos_archivo = encabezado + datos_ruido

    nombre_archivo = "archivo_falso.zip"
    guardar_archivo(nombre_archivo, datos_archivo)

    print(f"Archivo ZIP falso generado: {nombre_archivo} (tamaño: {len(datos_archivo)} bytes)")

if __name__ == "__main__":
    main()
