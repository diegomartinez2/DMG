import Markov_1
def generar_encabezado_zip():
    """
    Devuelve un encabezado típico de un archivo ZIP simple (puede ser ajustado).
    Los primeros 4 bytes suelen ser 'PK\x03\x04' indicando un archivo ZIP.
    """
    return bytes([0x50, 0x4B, 0x03, 0x04])

def main():
    estados = list(range(16))
    probabilidad_transicion = {
        0: {0: 0.7, 1: 0.1, 2: 0.1, 3: 0.1},
        1: {0: 0.6, 1: 0.2, 2: 0.1, 3: 0.1},
        2: {2: 0.8, 3: 0.1, 4: 0.1},
        3: {3: 0.7, 4: 0.3},
    }
    for estado in estados:
        if estado not in probabilidad_transicion:
            probabilidad_transicion[estado] = {e: 1/len(estados) for e in estados}

    matriz = crear_matriz_transicion(estados, probabilidad_transicion)
    estado_inicial = random.choice(estados)

    encabezado = generar_encabezado_zip()
    tamaño_total = 154
    longitud_ruido = tamaño_total - len(encabezado)

    datos_ruido = generar_cadena_markov(matriz, estado_inicial, longitud_ruido)

    datos_archivo = encabezado + datos_ruido

    nombre_archivo = "archivo_falso.zip"
    guardar_archivo(nombre_archivo, datos_archivo)

    print(f"Archivo ZIP simulado generado: {nombre_archivo}, tamaño: {len(datos_archivo)} bytes")

if __name__ == "__main__":
    main()
