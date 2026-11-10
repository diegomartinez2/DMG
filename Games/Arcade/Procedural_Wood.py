#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Procedural_Wood.py
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
#  Generador procedural de textura de madera al estilo POV-Ray con NumPy y Arcade.
#
'''
Generación de una textura de madera usando un método procedural.
Usando Arcade y PIL.
'''
import numpy as np
import arcade
from PIL import Image

SCREEN_WIDTH = 600
SCREEN_HEIGHT = 600
SCREEN_TITLE = "Textura Procedural de Madera - Estilo POV-Ray"

class WoodTextureWindow(arcade.Window):
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        arcade.set_background_color(arcade.color.BLACK)

        # 1. Crear el sprite que ocupará toda la pantalla
        self.texture_sprite = arcade.Sprite()
        self.texture_sprite.center_x = SCREEN_WIDTH / 2
        self.texture_sprite.center_y = SCREEN_HEIGHT / 2
        self.texture_sprite.width = SCREEN_WIDTH
        self.texture_sprite.height = SCREEN_HEIGHT

        # 2. Generar la textura de madera UNA SOLA VEZ al iniciar (adiós al lastre por frame)
        wood_image = self.generate_wood_texture(SCREEN_WIDTH, SCREEN_HEIGHT)

        # 3. Asignar la textura al sprite
        dynamic_texture = arcade.Texture(image=wood_image)
        self.texture_sprite.texture = dynamic_texture

        # 4. SpriteList para dibujar en Arcade 3.x
        self.sprite_list = arcade.SpriteList()
        self.sprite_list.append(self.texture_sprite)

    def generate_wood_texture(self, width, height):
        print("Calculando vetas de madera mediante matemáticas de vieja escuela...")

        # Crear malla de coordenadas centrada en la pantalla
        resolution_x = width
        resolution_y = height
        x_vals = np.linspace(-10, 10, resolution_x)
        y_vals = np.linspace(-10, 10, resolution_y)
        x, y = np.meshgrid(x_vals, y_vals)

        # TRUCO 2: Turbulencia matemática (deformar el espacio con senos cruzados)
        # Esto evita que los anillos sean círculos concéntricos aburridos de geometría pura
        turbulence = 0.8 * np.sin(x * 0.5) + 0.8 * np.cos(y * 0.5)

        # TRUCO 1: Coordenadas radiales estiradas para simular el tronco / tabla longitudinal
        # Estiramos el eje Y para darle aspecto de tabla alargada en vez de tronco redondo
        x_warp = x + turbulence
        y_warp = (y * 0.3) + turbulence
        r = np.sqrt(x_warp**2 + y_warp**2)

        # TRUCO 3: Frecuencia de los anillos y perfil de veta afilada
        # Usamos valor absoluto sobre un seno para crear crestas afiladas (vetas oscuras)
        frequency = 4.0
        pattern = np.abs(np.sin(r * frequency))

        # Añadir un armónico secundario para dar microdetalles en la madera
        pattern = 0.7 * pattern + 0.3 * np.abs(np.sin(r * frequency * 2.5 + x * 0.2))

        # Normalizar el patrón estrictamente entre 0.0 y 1.0
        p_min, p_max = pattern.min(), pattern.max()
        normalized = (pattern - p_min) / (p_max - p_min)

        # TRUCO 4: Paleta de colores RGB (Mapeo de tonos de madera)
        # Definimos color oscuro (vetas) y color claro (fondo de madera)
        # Color oscuro: Marrón café (Ej: RGB 40, 20, 10)
        # Color claro:  Marrón miel / caoba (Ej: RGB 180, 115, 60)
        color_dark = np.array([45, 24, 10], dtype=float)
        color_light = np.array([170, 110, 55], dtype=float)

        # Interpolar linealmente entre el color oscuro y el claro usando el patrón matemático
        # Expandimos 'normalized' a 3 dimensiones para multiplicar por RGB: (H, W, 1)
        factor = normalized[..., np.newaxis]
        rgb_matrix = color_dark * (1.0 - factor) + color_light * factor

        # Asegurar tipos de datos correctos para Pillow (uint8)
        img_array = np.clip(rgb_matrix, 0, 255).astype(np.uint8)

        # Convertir la matriz numpy a imagen PIL en formato RGBA
        pil_image = Image.fromarray(img_array, mode="RGB").convert("RGBA")
        return pil_image

    def on_update(self, delta_time: float):
        # Como es estática, aquí no hay que recalcular nada en cada frame. ¡Rendimiento puro!
        pass

    def on_draw(self):
        self.clear()
        self.sprite_list.draw()

def main():
    window = WoodTextureWindow()
    arcade.run()

if __name__ == "__main__":
    main()
