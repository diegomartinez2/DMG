#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Isometric_Grid.py
#  Renderizado básico de una retícula isométrica con Python y Arcade.
#

import arcade

SCREEN_WIDTH = 800
SCREEN_HEIGHT = 600
SCREEN_TITLE = "Proyección Isométrica - Estilo Clásico"

class IsometricWindow(arcade.Window):
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        self.background_color = arcade.color.DARK_SLATE_GRAY

        # Dimensiones de los tiles isométricos (proporción clásica 2:1)
        self.tile_width = 64
        self.tile_height = 32

        # Desplazamiento de la cámara para centrar el mapa en pantalla
        self.offset_x = SCREEN_WIDTH // 2
        self.offset_y = 150

    def grid_to_iso(self, gx, gy):
        """Convierte coordenadas de matriz lógica a coordenadas de pantalla isométrica."""
        screen_x = self.offset_x + (gx - gy) * (self.tile_width // 2)
        screen_y = self.offset_y + (gx + gy) * (self.tile_height // 2)
        return screen_x, screen_y

    def on_draw(self):
        self.clear()

        # Dibujar una cuadrícula isométrica de 8x8 tiles
        grid_size = 8
        for gx in range(grid_size):
            for gy in range(grid_size):
                # Obtener la posición isométrica del tile
                x, y = self.grid_to_iso(gx, gy)

                # Dibujar los rombos del suelo (simulando los tiles)
                # Como Arcade no dibuja polígonos arbitrarios complejos en una sola línea fácil sin ShapeList,
                # dibujamos líneas que conectan las esquinas del rombo o un círculo de referencia.

                # Coordenadas de los 4 vértices del rombo isométrico
                half_w = self.tile_width // 2
                half_h = self.tile_height // 2

                points = [
                    (x, y + half_h),          # Norte
                    (x + half_w, y),          # Este
                    (x, y - half_h),          # Sur
                    (x - half_w, y)           # Oeste
                ]

                # Dibujar el contorno del tile isométrico
                arcade.draw_polygon_outline(points, arcade.color.LIGHT_GRAY, 1)

                # Opcional: Marcar el centro del tile con un punto
                arcade.draw_circle_filled(x, y, 3, arcade.color.CRIMSON)

def main():
    window = IsometricWindow()
    arcade.run()

if __name__ == "__main__":
    main()
