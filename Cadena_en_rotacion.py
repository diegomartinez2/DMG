import numpy as np
import matplotlib.pyplot as plt

def calcular_forma_cadena_rotacion(
    num_puntos,
    distancia_eje,
    factor_rigidez,
    num_iteraciones
):
    """
    Calcula la forma de una cadena "colgada" por la fuerza centrífuga.

    La cadena está sujeta a dos puntos en el eje de rotación. La forma de la
    cadena se determina por la fuerza centrífuga, que aumenta con la distancia
    al eje de rotación.

    Parámetros:
    num_puntos (int): Número de puntos discretos que representan la cadena.
    distancia_eje (float): Distancia entre los dos puntos de sujeción a lo largo del eje.
    factor_rigidez (float): Un parámetro que representa la relación entre la fuerza
                             centrífuga y la tensión de la cadena. Un valor más alto
                             significa una curva más pronunciada.
    num_iteraciones (int): Número de pasos de relajación para que la curva converja.

    Retorna:
    tuple: Una tupla con las coordenadas x e y de la curva.
    """

    # Coordenadas y: se mantienen fijas a lo largo del eje de rotación
    y_coords = np.linspace(-distancia_eje / 2, distancia_eje / 2, num_puntos)

    # Coordenadas x: se inicializan a cero (cadena recta) y se actualizan
    x_coords = np.zeros(num_puntos)

    # Las coordenadas de los extremos de la cadena se fijan en 0, ya que están
    # en el eje de rotación
    x_coords[0] = 0
    x_coords[-1] = 0

    # Iteramos para encontrar la forma de equilibrio
    # El método de relajación actualiza la posición de cada punto en función de
    # sus vecinos y la fuerza centrífuga
    for _ in range(num_iteraciones):
        for i in range(1, num_puntos - 1):
            # Ecuación de diferencias finitas simplificada para la fuerza centrífuga
            # La fórmula es: x_new = (x_prev + x_next) / (2 + C)
            # Donde C es el factor de rigidez que depende de la tensión, masa y velocidad angular
            x_coords[i] = (x_coords[i-1] + x_coords[i+1]) / (2 + factor_rigidez)

    return x_coords, y_coords

# --- Configuración de la simulación ---
num_puntos = 100
distancia_eje = 10.0 # Metros
factor_rigidez = 0.5 # Ajusta este valor para ver diferentes curvas
num_iteraciones = 5000

# Calcular la forma de la curva
x, y = calcular_forma_cadena_rotacion(num_puntos, distancia_eje, factor_rigidez, num_iteraciones)

# --- Visualización del resultado ---
plt.figure(figsize=(8, 6))
plt.plot(x, y, label='Forma de la cadena')
plt.xlabel('Distancia al eje de rotación (m)')
plt.ylabel('Posición a lo largo del eje (m)')
plt.title('Forma de una cadena suspendida por rotación')
plt.grid(True)
plt.axvline(0, color='gray', linestyle='--', linewidth=0.8) # Dibuja el eje de rotación
plt.axhline(y[0], color='gray', linestyle='--', linewidth=0.8) # Puntos de sujeción
plt.axhline(y[-1], color='gray', linestyle='--', linewidth=0.8)
plt.legend()
plt.show()
