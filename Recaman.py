def secuencia_recaman(n_terminos):
    """Genera los primeros n términos de la secuencia de Recamán.

    :param n_terminos: Cantidad de números de la secuencia a calcular.
    :return: Lista con la secuencia de Recamán.
    """
    if n_terminos <= 0:
        return []

    # Inicializamos la secuencia con el primer término (A0 = 0)
    secuencia = [0]

    # Usamos un conjunto (set) para saber qué números ya han salido.
    # Buscar en un set en Python toma tiempo O(1), es muchísimo más rápido que buscar en una lista.
    visitados = {0}

    # El bucle empieza desde el paso 1 hasta el paso n_terminos - 1
    for paso in range(1, n_terminos):
        # Obtenemos el último número que añadimos a la secuencia
        ultimo_numero = secuencia[-1]

        # REGLA 1: Intentamos dar un paso hacia atrás (restar el paso actual)
        atras = ultimo_numero - paso

        # REGLA 2: Evaluamos si podemos quedarnos con el número de atrás.
        # Debe ser estrictamente mayor que cero Y no haber salido antes.
        if atras > 0 and atras not in visitados:
            proximo_numero = atras
        else:
            # Si no cumple, obligatoriamente damos un paso hacia adelante (sumar el paso)
            proximo_numero = ultimo_numero + paso

        # Guardamos el resultado en la lista y en el conjunto de control
        secuencia.append(proximo_numero)
        visitados.add(proximo_numero)

    return secuencia


# --- PRUEBA DEL SCRIPT ---
if __name__ == "__main__":
    # Calculamos los primeros 15 términos como ejemplo
    cantidad = 15
    resultado = secuencia_recaman(cantidad)

    print(f"--- Los primeros {cantidad} términos de Recamán ---")
    for i, valor in enumerate(resultado):
        print(f"A({i}) = {valor}")
