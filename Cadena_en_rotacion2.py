import numpy as np
import matplotlib.pyplot as plt

def calcular_curva_centrifuga(L, N, omega, rho, num_iteraciones, dt):
    """
    Calcula la forma de una curva suspendida en rotación usando el método de relajación.

    Parámetros:
    L (float): Longitud total de la cadena.
    N (int): Número de puntos discretos en la cadena (incluyendo los extremos).
    omega (float): Velocidad angular de rotación (rad/s).
    rho (float): Densidad de masa lineal de la cadena (kg/m).
    num_iteraciones (int): Número de iteraciones para la relajación.
    dt (float): Factor de "paso de tiempo" o de amortiguación para las actualizaciones.

    Retorna:
    tuple: Arrays de las coordenadas x e y de la curva.
    """
    # Inicialización de los puntos.
    # Los puntos de anclaje están en el eje de rotación (x=0).
    # La cadena se inicia como una línea recta.
    x = np.zeros(N)
    y = np.linspace(-L / 2, L / 2, N)

    # La masa de cada punto.
    masa_punto = (L / (N - 1)) * rho

    # Magnitud de la tensión, es un valor que debe ser sintonizado o calculado
    # más rigurosamente. Para este ejemplo, lo asumimos como un valor fijo.
    # Cuanto mayor sea la tensión, más "rígida" es la cadena.
    tension_magnitud = 100.0

    # Bucle principal de relajación
    for _ in range(num_iteraciones):
        x_new = np.copy(x)
        y_new = np.copy(y)

        # Iteramos sobre los puntos internos (excluyendo los extremos fijos)
        for i in range(1, N - 1):
            # Vector del punto anterior al actual
            vec_prev = np.array([x[i-1] - x[i], y[i-1] - y[i]])
            # Vector del punto siguiente al actual
            vec_next = np.array([x[i+1] - x[i], y[i+1] - y[i]])

            # Fuerzas de tensión de los segmentos adyacentes
            fuerza_tension = (tension_magnitud * (vec_prev / np.linalg.norm(vec_prev)) +
                             tension_magnitud * (vec_next / np.linalg.norm(vec_next)))

            # Fuerza centrípeta
            fuerza_centripeta = np.array([masa_punto * omega**2 * x[i], 0.0])

            # Fuerza neta sobre el punto
            fuerza_neta = fuerza_tension + fuerza_centripeta

            # Actualización de la posición del punto
            x_new[i] += dt * fuerza_neta[0]
            y_new[i] += dt * fuerza_neta[1]

        # Actualizamos las posiciones de todos los puntos
        x = x_new
        y = y_new

    return x, y

if __name__ == "__main__":
    # Parámetros del sistema
    longitud_cadena = 10.0  # metros
    num_puntos = 101
    velocidad_angular = 2.0  # rad/s
    densidad_masa_lineal = 1.0  # kg/m

    # Parámetros del método numérico
    iteraciones = 1000
    paso_tiempo = 0.05

    # Calcular la curva
    x_curva, y_curva = calcular_curva_centrifuga(
        longitud_cadena, num_puntos, velocidad_angular, densidad_masa_lineal, iteraciones, paso_tiempo
    )

    # Graficar la curva
    plt.figure(figsize=(8, 6))
    plt.plot(x_curva, y_curva, 'o-', markersize=2, label='Cadena en rotación')
    plt.plot([x_curva[0], x_curva[-1]], [y_curva[0], y_curva[-1]], 'ro', label='Puntos de anclaje')
    plt.title('Forma de la Cadena bajo Fuerza Centrípeta')
    plt.xlabel('Distancia al Eje de Rotación (x) [m]')
    plt.ylabel('Posición a lo largo del Eje de Rotación (y) [m]')
    plt.grid(True)
    plt.legend()
    plt.axis('equal')  # Mantiene la proporción para que la curva no se vea distorsionada
    plt.show()
