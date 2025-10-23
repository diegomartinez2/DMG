import numpy as np
import matplotlib.pyplot as plt

def mantener_longitud_segmentos(x, y, L_segmento):
    """
    Ajusta las posiciones de los puntos (excepto los extremos fijos)
    para mantener la longitud de cada segmento.
    """
    N = len(x)

    # 1. Corrección hacia adelante (desde el anclaje inicial x[0], y[0])
    # Este bucle ajusta i basándose en la posición de i-1
    for i in range(1, N - 1):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        distancia_actual = np.sqrt(dx**2 + dy**2)

        # Evitar división por cero si los puntos están en el mismo lugar
        if distancia_actual == 0:
            continue

        factor_correccion = L_segmento / distancia_actual
        # Mueve x[i], y[i] para que esté a L_segmento de x[i-1], y[i-1]
        x[i] = x[i-1] + dx * factor_correccion
        y[i] = y[i-1] + dy * factor_correccion

    # 2. Corrección hacia atrás (desde el anclaje final x[N-1], y[N-1])
    # Este bucle ajusta i basándose en la posición de i+1
    # Esto es crucial para que la restricción se propague desde ambos extremos
    for i in range(N - 2, 0, -1):
        dx = x[i] - x[i+1]
        dy = y[i] - y[i+1]
        distancia_actual = np.sqrt(dx**2 + dy**2)

        if distancia_actual == 0:
            continue

        factor_correccion = L_segmento / distancia_actual
        # Mueve x[i], y[i] para que esté a L_segmento de x[i+1], y[i+1]
        x[i] = x[i+1] + dx * factor_correccion
        y[i] = y[i+1] + dy * factor_correccion

    return x, y

def calcular_curva_centrifuga_con_correccion(L_anclajes, N, omega, rho, num_iteraciones, dt):
    """
    Calcula la forma de la curva usando un método de relajación.
    L_anclajes es la distancia en 'y' entre los puntos de anclaje.
    """

    # --- 1. Inicialización ---
    # Inicialización de los puntos con una forma parabólica
    # Los anclajes están en el eje Y (x=0) en y = -L/2 y y = +L/2
    y = np.linspace(-L_anclajes / 2, L_anclajes / 2, N)

    sag_inicial = 2.0  # "Comba" inicial de la cadena
    a = -4 * sag_inicial / (L_anclajes**2)
    # x=0 en los extremos, x=sag_inicial en el centro
    x = a * y**2 + sag_inicial

    # Guardar posiciones de anclaje (que nunca deben moverse)
    x_anclaje_0, y_anclaje_0 = x[0], y[0]
    x_anclaje_N, y_anclaje_N = x[-1], y[-1]

    # --- 2. Calcular propiedades de la cadena ---
    # Calcular la longitud total real de la curva parabólica inicial
    L_total_cadena = 0
    for i in range(N - 1):
        L_total_cadena += np.sqrt((x[i+1] - x[i])**2 + (y[i+1] - y[i])**2)

    # Esta es la longitud de segmento que DEBE mantenerse
    L_segmento = L_total_cadena / (N - 1)

    masa_punto = L_segmento * rho

    # --- 3. Simulación (Método de Relajación) ---
    for _ in range(num_iteraciones):

        # Copiamos el estado actual para no interferir con los cálculos
        x_new = np.copy(x)
        y_new = np.copy(y)

        # --- A. Aplicar fuerzas EXTERNAS (Centrífuga) ---
        # Solo iteramos sobre los puntos interiores (1 a N-2)
        for i in range(1, N - 1):
            # Fuerza centrífuga: F = m * omega^2 * r
            # El radio 'r' es la distancia al eje Y, que es x[i]
            # La fuerza es en la dirección +x
            f_centrifuga_x = masa_punto * (omega**2) * x[i]
            f_centrifuga_y = 0.0

            # Aplicamos un "empujón" en la dirección de la fuerza
            # 'dt' aquí no es un paso de tiempo real, sino un
            # factor de "relajación" o "aprendizaje".
            x_new[i] += dt * f_centrifuga_x
            y_new[i] += dt * f_centrifuga_y

        # Actualizamos las posiciones con el "empujón"
        x, y = x_new, y_new

        # --- B. Aplicar Restricciones (Tensión Implícita) ---
        # Forzamos a los puntos a volver a sus longitudes de segmento correctas.
        # Esto simula la tensión tirando de los puntos.
        x, y = mantener_longitud_segmentos(x, y, L_segmento)

        # (Opcional pero seguro) Re-forzar las posiciones de los anclajes
        # por si la función de restricción tuviera algún error.
        x[0], y[0] = x_anclaje_0, y_anclaje_0
        x[-1], y[-1] = x_anclaje_N, y_anclaje_N


    return x, y

if __name__ == "__main__":
    # Parámetros
    longitud_eje_anclajes = 10.0 # Distancia 'L' entre anclajes en el eje Y
    num_puntos = 101             # Número de puntos (N) para discretizar la cadena
    velocidad_angular = 2.0      # omega (rad/s)
    densidad_masa_lineal = 1.0   # rho (kg/m)

    iteraciones = 2000           # Aumentado para mejor convergencia
    paso_tiempo_relax = 0.005    # 'dt' (factor de relajación, ajustado para estabilidad)

    print(f"Iniciando simulación...")

    x_curva, y_curva = calcular_curva_centrifuga_con_correccion(
        longitud_eje_anclajes, num_puntos, velocidad_angular, densidad_masa_lineal,
        iteraciones, paso_tiempo_relax
    )

    print("Simulación completada. Mostrando gráfico.")

    # --- Graficar ---
    plt.figure(figsize=(8, 8))
    plt.plot(x_curva, y_curva, 'o-', markersize=2, label=f'Cadena en rotación ($\\omega$={velocidad_angular} rad/s)')
    plt.plot([x_curva[0], x_curva[-1]], [y_curva[0], y_curva[-1]], 'ro', markersize=8, label='Puntos de anclaje')
    plt.title('Forma de Equilibrio de Cadena Rotatoria (Sin Gravedad)')
    plt.xlabel('Distancia al Eje de Rotación (x) [m]')
    plt.ylabel('Posición a lo largo del Eje de Rotación (y) [m]')
    plt.grid(True)
    plt.legend()
    # plt.axis('equal') es crucial para ver la forma real
    plt.axis('equal')
    plt.show()
