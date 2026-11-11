#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Pendulo.py
#  Simulación física de un péndulo mediante matemáticas puras (Ecuación diferencial).
#

import math
import arcade

SCREEN_WIDTH = 600
SCREEN_HEIGHT = 600
SCREEN_TITLE = "Simulación de Péndulo - Matemáticas Puras"

class PendulumWindow(arcade.Window):
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        self.background_color = arcade.color.DARK_SLATE_GRAY

        # Punto de pivote fijo en la parte superior central
        self.pivot_x = SCREEN_WIDTH // 2
        self.pivot_y = SCREEN_HEIGHT - 100

        # Variables físicas del péndulo
        self.length = 250.0          # Longitud de la cuerda en píxeles
        self.angle = math.pi / 3     # Ángulo inicial (60 grados en radianes)
        self.angular_velocity = 0.0  # Velocidad angular inicial

        # Parámetros de simulación
        self.gravity = 1500.0        # Gravedad simulada adaptada a escala de píxeles
        self.damping = 0.998         # Factor de fricción del aire (amortiguamiento)

    def on_update(self, delta_time: float):
        # 1. Calcular la aceleración angular usando la ecuación diferencial del péndulo simple:
        #    alpha = -(g / L) * sin(theta)
        angular_acceleration = -(self.gravity / self.length) * math.sin(self.angle)

        # 2. Integración numérica (Método de Euler semi-implícito)
        self.angular_velocity += angular_acceleration * delta_time
        self.angular_velocity *= self.damping  # Aplicar pérdida de energía por fricción
        self.angle += self.angular_velocity * delta_time

    def on_draw(self):
        # Limpiar la pantalla
        self.clear()

        # 3. Calcular las coordenadas cartesianas (x, y) de la masa usando trigonometría pura
        #    El ángulo 0 es totalmente vertical hacia abajo.
        x = self.pivot_x + self.length * math.sin(self.angle)
        y = self.pivot_y - self.length * math.cos(self.angle)

        # Dibujar la cuerda del péndulo
        arcade.draw_line(
            self.pivot_x, self.pivot_y,
            x, y,
            arcade.color.LIGHT_GRAY,
            line_width=3
        )

        # Dibujar el perno de sujeción superior
        arcade.draw_circle_filled(self.pivot_x, self.pivot_y, 8, arcade.color.WHITE)

        # Dibujar la masa oscilante en el extremo
        arcade.draw_circle_filled(x, y, 22, arcade.color.CRIMSON)
        arcade.draw_circle_outline(x, y, 22, arcade.color.WHITE, 2)

def main():
    window = PendulumWindow()
    arcade.run()

if __name__ == "__main__":
    main()
