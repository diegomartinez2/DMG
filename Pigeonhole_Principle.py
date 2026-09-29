#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Pingeon_Principle.py
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
from collections import deque

def encontrar_factor_binario(N: int):
    """
    Dado un entero positivo N, encuentra el factor k tal que (N * k)
    sea el menor número entero positivo formado únicamente por los dígitos '0' y '1'.
    """
    if N <= 0:
        raise ValueError("El número debe ser un entero positivo mayor que 0.")

    if N == 1:
        return 1, 1

    # parent[resto] almacenará (resto_anterior, dígito_usado)
    # Permite reconstruir el número final sin guardar cadenas gigantes en memoria.
    parent = {}
    queue = deque()

    # Empezamos la búsqueda construyendo el número desde el primer dígito '1'
    start_rem = 1 % N
    parent[start_rem] = (-1, '1')
    queue.append(start_rem)

    encontrado = False

    while queue:
        curr_rem = queue.popleft()

        if curr_rem == 0:
            encontrado = True
            break

        # Transiciones agregando '0' o '1' al final del número actual
        for digit in (0, 1):
            next_rem = (curr_rem * 10 + digit) % N

            # Si es la primera vez que visitamos este resto, lo registramos
            if next_rem not in parent:
                parent[next_rem] = (curr_rem, str(digit))
                queue.append(next_rem)

                if next_rem == 0:
                    encontrado = True
                    break
        if encontrado:
            break

    # Reconstrucción del número de 0s y 1s recorriendo hacia atrás los restos
    digits = []
    curr = 0
    while curr != -1:
        prev_rem, d = parent[curr]
        digits.append(d)
        curr = prev_rem

    digits.reverse()
    numero_binario_str = "".join(digits)
    numero_binario = int(numero_binario_str)
    factor = numero_binario // N

    return factor, numero_binario


if __name__ == "__main__":
    # Prueba con múltiplos de 9 y números complejos
    casos = [9, 99, 999, 7, 13, 90]

    for n in casos:
        k, res = encontrar_factor_binario(n)
        print(f"N = {n}")
        print(f"Factor k = {k}")
        print(f"N * k = {res}")
        print("-" * 50)
