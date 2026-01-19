#!/usr/bin/env python

import hashlib
import itertools
import string
import time
from typing import Optional, Set

# --- Parámetros de Configuración ---

# El hash objetivo (debe ser generado a partir de una palabra conocida para la prueba)
# Ejemplo: 'password' -> c5a6042c13d94066914619736c072c418c3445e4125f16212ac2d951a84f506e
HASH_OBJETIVO = "c5a6042c13d94066914619736c072c418c3445e4125f16212ac2d951a84f506e"

# Conjunto de caracteres a probar.
# ADVERTENCIA: Aumentar este conjunto o la longitud aumenta exponencialmente el tiempo.
CARACTERES_A_PROBAR = string.ascii_lowercase + string.digits # a-z y 0-9

# Rango de longitudes de palabra a probar
LONGITUD_MIN = 1
LONGITUD_MAX = 7 # ¡NO aumentar mucho más allá de 7! Longitudes mayores tardarán horas/días/años.

# --- Funciones de Utilidad ---

def generar_hash_sha256(palabra: str) -> str:
    """Calcula el hash SHA-256 de una palabra y lo retorna en formato hexadecimal."""
    return hashlib.sha256(palabra.encode('utf-8')).hexdigest()

def ataque_fuerza_bruta(hash_objetivo: str, min_len: int, max_len: int, chars: str) -> Optional[str]:
    """
    Intenta encontrar la palabra original probando todas las combinaciones posibles.
    """
    print(f"Iniciando ataque de Fuerza Bruta contra el hash: {hash_objetivo[:16]}...")
    print(f"Probando caracteres: {chars[:20]}...")

    start_time = time.time()
    total_intentos = 0

    for length in range(min_len, max_len + 1):
        print(f"\n--- Probando palabras de longitud {length} ---")

        # Generar todas las combinaciones posibles de esta longitud
        for intento_tuple in itertools.product(chars, repeat=length):
            palabra_intento = "".join(intento_tuple)
            total_intentos += 1

            # Calcular el hash
            hash_intento = generar_hash_sha256(palabra_intento)

            # Comprobar la coincidencia
            if hash_intento == hash_objetivo:
                end_time = time.time()
                tiempo_transcurrido = end_time - start_time
                print("\n" + "="*50)
                print(f"¡ÉXITO! Palabra encontrada:")
                print(f"Palabra: {palabra_intento}")
                print(f"Hash: {hash_intento}")
                print(f"Tiempo Total: {tiempo_transcurrido:.2f} segundos")
                print(f"Intentos: {total_intentos:,}")
                print("="*50)
                return palabra_intento

            # Mostrar progreso cada cierto número de intentos
            if total_intentos % 1_000_000 == 0:
                elapsed = time.time() - start_time
                print(f"  > Intentos: {total_intentos:,} | Velocidad: {total_intentos/elapsed:.2f} hashes/s", end='\r')

    # Si se termina el bucle sin éxito
    end_time = time.time()
    print("\n" + "="*50)
    print("FALLO: La palabra no se encontró dentro de los parámetros definidos.")
    print(f"Intentos Totales: {total_intentos:,}")
    print(f"Tiempo Total: {end_time - start_time:.2f} segundos")
    print("Considere aumentar HASH_OBJETIVO (o el conjunto de caracteres), pero el tiempo aumentará exponencialmente.")
    print("="*50)
    return None

# --- Ejecución Principal ---

if __name__ == "__main__":
    # La palabra que corresponde al HASH_OBJETIVO es 'password'
    # Ejecuta el ataque
    ataque_fuerza_bruta(
        hash_objetivo=HASH_OBJETIVO,
        min_len=LONGITUD_MIN,
        max_len=LONGITUD_MAX,
        chars=CARACTERES_A_PROBAR
    )
