def fibonacci_python(n):
    """Genera la secuencia de Fibonacci hasta el n-ésimo término."""
    secuencia = []
    a, b = 0, 1

    for _ in range(n):
        secuencia.append(a)
        # Asignación múltiple (característica común con Julia)
        a, b = b, a + b

    return secuencia

# Ejecución
n_terminos = 10
resultado = fibonacci_python(n_terminos)
print(f"Python (n={n_terminos}): {resultado}")
