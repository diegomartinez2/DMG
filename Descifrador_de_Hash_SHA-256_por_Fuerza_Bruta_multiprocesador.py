#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Descifrador_de_Hash_SHA-256_por_Fuerza_Bruta.py
#
#  Copyright 2026 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
# ---------------------------

import hashlib
import itertools
import multiprocessing
import string
import time
from typing import Optional, Tuple

# --- Parámetros de Configuración ---

# El hash objetivo (debe ser generado a partir de una palabra conocida para la prueba)
# Ejemplo: 'password' -> c5a6042c13d94066914619736c072c418c3445e4125f16212ac2d951a84f506e
HASH_OBJETIVO = "c5a6042c13d94066914619736c072c418c3445e4125f16212ac2d951a84f506e"
CARACTERES_A_PROBAR = string.ascii_lowercase + string.digits  # a-z y 0-9
LONGITUD_MIN = 1
LONGITUD_MAX = 7
NUM_PROCESOS = 4  # Número de núcleos de CPU a utilizar


def generar_hash_sha256(palabra: str) -> str:
    """Calcula el hash SHA-256 de una cadena en hexadecimal."""
    return hashlib.sha256(palabra.encode('utf-8')).hexdigest()


def probar_rango_prefijo(args: Tuple[str, str, int, int]) -> Optional[str]:
    """
    Función que ejecuta cada proceso trabajador.
    Prueba combinaciones que comienzan con un prefijo específico.
    """
    prefijo, chars, min_len, max_len = args

    for length in range(min_len, max_len + 1):
        # Si la longitud es 1, solo se prueba el propio prefijo
        if length == 1:
            if generar_hash_sha256(prefijo) == HASH_OBJETIVO:
                return prefijo
            continue

        # Genera el resto de los caracteres de la combinación
        resto_longitud = length - 1
        for intento_tuple in itertools.product(chars, repeat=resto_longitud):
            palabra_intento = prefijo + "".join(intento_tuple)
            if generar_hash_sha256(palabra_intento) == HASH_OBJETIVO:
                return palabra_intento

    return None


def ataque_fuerza_bruta_multiproceso() -> Optional[str]:
    """Coordina los 4 procesos para realizar la búsqueda en paralelo."""
    print(f"Iniciando búsqueda paralela utilizando {NUM_PROCESOS} núcleos de CPU...")
    start_time = time.time()

    # Se crean las tareas asignando cada carácter inicial como trabajo independiente
    tareas = [
        (char_inicial, CARACTERES_A_PROBAR, LONGITUD_MIN, LONGITUD_MAX)
        for char_inicial in CARACTERES_A_PROBAR
    ]

    # Crear el pool de procesos (limitado a 4 núcleos)
    with multiprocessing.Pool(processes=NUM_PROCESOS) as pool:
        # imap_unordered procesa los resultados conforme los trabajadores van terminando
        for resultado in pool.imap_unordered(probar_rango_prefijo, tareas):
            if resultado is not None:
                # Si un proceso encuentra la coincidencia, cancela los demás trabajadores
                pool.terminate()
                end_time = time.time()

                print("\n" + "=" * 50)
                print(f"¡ÉXITO! Palabra encontrada:")
                print(f"Palabra: {resultado}")
                print(f"Hash: {generar_hash_sha256(resultado)}")
                print(f"Tiempo Total: {end_time - start_time:.2f} segundos")
                print("=" * 50)
                return resultado

    end_time = time.time()
    print("\n" + "=" * 50)
    print("FALLO: La palabra no se encontró dentro de los parámetros definidos.")
    print(f"Tiempo Total: {end_time - start_time:.2f} segundos")
    print("=" * 50)
    return None


if __name__ == "__main__":
    ataque_fuerza_bruta_multiproceso()
