#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Double_Pendulum.py
#  Simulación de Péndulo Doble caótico con Mecánica Lagragiana y Runge-Kutta 4 en Arcade.
#

import math
import arcade
from collections import deque

SCREEN_WIDTH = 800
SCREEN_HEIGHT = 800
SCREEN_TITLE = "Simulación de Péndulo Doble - Caos Determinista"

class DoublePendulumWindow(arcade.Window):
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        self.background_color = arcade.color.DARK_SLATE_GRAY

        # Pivote en el centro superior de la pantalla
        self.pivot_x = SCREEN_WIDTH // 2
        self.pivot_y = SCREEN_HEIGHT - 200

        # Parámetros físicos
        self.L1 = 180.0  # Longitud de la barra 1
        self.L2 = 180.0  # Longitud de la barra 2
        self.m1 = 10.0   # Masa 1
        self.m2 = 10.0   # Masa 2
        self.g = 981.0   # Gravedad simulada

        # Estado inicial: [theta1, omega1, theta2, omega2]
        # Empezamos con ángulos elevados para buscar el caos inmediato
        self.state = [math.pi / 2, 0.0, math.pi / 2, 0.0]

        # Estela (historial de posiciones de la segunda masa para pintar el caos)
        self.trail = deque(maxlen=300)

    def derivatives(self, state):
        """
        Ecuaciones de movimiento derivadas de Lagrange para el péndulo doble.
        Devuelve las derivadas de las variables de estado: [d_theta1, d_omega1, d_theta2, d_omega2]
        """
        t1, w1, t2, w2 = state
        g = self.g
        l1, l2 = self.L1, self.L2
        m1, m2 = self.m1, self.m2

        delta = t1 - t2

        # Denominador común de las ecuaciones de Lagrange
        den = l1 * (2 * m1 + m2 - m2 * math.cos(2 * t1 - 2 * t2))

        # Aceleración angular 1 (alpha 1)
        num1 = -g * (2 * m1 + m2) * math.sin(t1) - m2 * g * math.sin(t1 - 2 * t2) - 2 * math.sin(delta) * m2 * (w2**2 * l2 + w1**2 * l1 * math.cos(delta))
        alpha1 = num1 / den

        # Aceleración angular 2 (alpha 2)
        num2 = 2 * math.sin(delta) * (w1**2 * l1 * (m1 + m2) + g * (m1 + m2) * math.cos(t1) + w2**2 * l2 * m2 * math.cos(delta))
        alpha2 = num2 / (l2 * (2 * m1 + m2 - m2 * math.cos(2 * t1 - 2 * t2)))

        return [w1, alpha1, w2, alpha2]

    def rk4_step(self, dt):
        """
        Paso de integración numérica Runge-Kutta de 4º orden (RK4).
        Garantiza que la energía del péndulo no se descontrole.
        """
        y = self.state

        # k1
        k1 = self.derivatives(y)

        # k2
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(4)]
        k2 = self.derivatives(y2)

        # k3
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(4)]
        k3 = self.derivatives(y3)

        # k4
        y4 = [y[i] + dt * k3[i] for i in range(4)]
        k4 = self.derivatives(y4)

        # Actualizar estado sumando las 4 pendientes ponderadas
        for i in range(4):
            y[i] += (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])

    def on_update(self, delta_time: float):
        # Limitar el delta_time máximo para evitar saltos raros si la CPU se atraganta
        dt = min(delta_time, 0.05)

        # Aplicamos varios sub-pasos por frame para mayor precisión numérica
        substeps = 5
        for _ in range(substeps):
            self.rk4_step(dt / substeps)

        # Calcular coordenadas cartesianas actuales para el dibujo y la estela
        t1, _, t2, _ = self.state
        x1 = self.pivot_x + self.L1 * math.sin(t1)
        y1 = self.pivot_y - self.L1 * math.cos(t1)
        x2 = x1 + self.L2 * math.sin(t2)
        y2 = y1 - self.L2 * math.cos(t2)

        # Guardar posición en la estela
        self.trail.append((x2, y2))

    def on_draw(self):
        self.clear()

        # 1. Calcular posiciones cartesianas finales
        t1, _, t2, _ = self.state
        x1 = self.pivot_x + self.L1 * math.sin(t1)
        y1 = self.pivot_y - self.L1 * math.cos(t1)
        x2 = x1 + self.L2 * math.sin(t2)
        y2 = y1 - self.L2 * math.cos(t2)

        # 2. Dibujar la estela caótica de la segunda masa
        if len(self.trail) > 1:
            points = list(self.trail)
            for i in range(len(points) - 1):
                # Desvanecer la estela según la antigüedad del punto
                alpha = int(255 * (i / len(points)))
                color = (255, 100, 100, alpha)
                arcade.draw_line(points[i][0], points[i][1], points[i+1][0], points[i+1][1], color, line_width=2)

        # 3. Dibujar las barras del péndulo
        arcade.draw_line(self.pivot_x, self.pivot_y, x1, y1, arcade.color.LIGHT_GRAY, line_width=4)
        arcade.draw_line(x1, y1, x2, y2, arcade.color.LIGHT_GRAY, line_width=4)

        # 4. Dibujar pivote, primera masa y segunda masa
        arcade.draw_circle_filled(self.pivot_x, self.pivot_y, 8, arcade.color.WHITE)
        arcade.draw_circle_filled(x1, y1, 14, arcade.color.CYAN)
        arcade.draw_circle_filled(x2, y2, 14, arcade.color.CRIMSON)
        arcade.draw_circle_outline(x2, y2, 14, arcade.color.WHITE, 2)

def main():
    window = DoublePendulumWindow()
    arcade.run()

if __name__ == "__main__":
    main()
