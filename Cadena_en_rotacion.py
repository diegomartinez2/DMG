import numpy as np
import matplotlib.pyplot as plt

def mantener_longitud_segmentos(x, y, L_segmento):
    """
    Ajusta la posición de los puntos para mantener una distancia constante
    entre los eslabones, manteniendo ambos extremos fijos.
    """
    N = len(x)

    # 1. Corrección hacia adelante (desde el primer anclaje)
    for i in range(1, N - 1):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        distancia_actual = np.sqrt(dx**2 + dy**2)

        factor_correccion = L_segmento / distancia_actual

        x[i] = x[i-1] + dx * factor_correccion
        y[i] = y[i-1] + dy * factor_correccion

    # 2. Corrección hacia atrás (desde el último anclaje)
    for i in range(N - 2, 0, -1):
        dx = x[i] - x[i+1]
        dy = y[i] - y[i+1]
        distancia_actual = np.sqrt(dx**2 + dy**2)

        factor_correccion = L_segmento / distancia_actual

        x[i] = x[i+1] + dx * factor_correccion
        y[i] = y[i+1] + dy * factor_correccion

    return x, y

def calcular_curva_centrifuga_con_correccion(L, N, omega, rho, num_iteraciones, dt):
    """
    Calcula la forma de la curva con corrección de la longitud de los segmentos.
    """
    L_segmento = L / (N - 1)

    # Inicialización de los puntos.
    y = np.linspace(-L / 2, L / 2, N)
    #x = np.zeros(N)
    # Calcular los coeficientes de una parábola que pasa por los puntos de anclaje (0, +/-L/2)
    # y tiene un "sag" inicial (distancia máxima al eje y)
    sag_inicial = 0.5  # Puedes ajustar este valor
    a = -4 * sag_inicial / (L**2)
    x = a * y**2 + sag_inicial
    masa_punto = L_segmento * rho
    tension_magnitud = 100.0

    for _ in range(num_iteraciones):
        x_new = np.copy(x)
        y_new = np.copy(y)

        # 1. Calcular y aplicar fuerzas (movimiento)
        for i in range(1, N - 1):
            vec_prev = np.array([x[i-1] - x[i], y[i-1] - y[i]])
            vec_next = np.array([x[i+1] - x[i], y[i+1] - y[i]])

            fuerza_tension = (tension_magnitud * (vec_prev / np.linalg.norm(vec_prev)) +
                             tension_magnitud * (vec_next / np.linalg.norm(vec_next)))

            fuerza_centripeta = np.array([masa_punto * omega**2 * x[i], 0.0])
            fuerza_neta = fuerza_tension + fuerza_centripeta

            x_new[i] += dt * fuerza_neta[0]
            y_new[i] += dt * fuerza_neta[1]

        x, y = x_new, y_new

        # 2. Corregir la longitud de los segmentos manteniendo ambos extremos fijos
        x, y = mantener_longitud_segmentos(x, y, L_segmento)

    return x, y

if __name__ == "__main__":
    # Parámetros del sistema
    longitud_cadena = 10.0  # metros
    num_puntos = 101
    velocidad_angular = 2.0  # rad/s
    densidad_masa_lineal = 1.0  # kg/m

    # Parámetros del método numérico
    iteraciones = 1000
    paso_tiempo = 0.01

    # Calcular la curva
    x_curva, y_curva = calcular_curva_centrifuga_con_correccion(
        longitud_cadena, num_puntos, velocidad_angular, densidad_masa_lineal, iteraciones, paso_tiempo
    )

    # Graficar la curva
    plt.figure(figsize=(8, 6))
    plt.plot(x_curva, y_curva, 'o-', markersize=2, label='Cadena en rotación (longitud constante)')
    plt.plot([x_curva[0], x_curva[-1]], [y_curva[0], y_curva[-1]], 'ro', label='Puntos de anclaje')
    plt.title('Forma de la Cadena bajo Fuerza Centrípeta')
    plt.xlabel('Distancia al Eje de Rotación (x) [m]')
    plt.ylabel('Posición a lo largo del Eje de Rotación (y) [m]')
    plt.grid(True)
    plt.legend()
    plt.axis('equal')
    plt.show()
