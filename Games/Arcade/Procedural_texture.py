#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Procedural_texture.py
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
#
import numpy as np
import arcade
from PIL import Image

SCREEN_WIDTH = 600
SCREEN_HEIGHT = 600
SCREEN_TITLE = "Textura Dinámica con NumPy y Arcade"

class ProceduralTextureWindow(arcade.Window):
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        arcade.set_background_color(arcade.color.BLACK)

        self.time_offset = 0.0

        # 1. Crear el sprite con las dimensiones de la pantalla
        self.texture_sprite = arcade.Sprite()
        self.texture_sprite.center_x = SCREEN_WIDTH / 2
        self.texture_sprite.center_y = SCREEN_HEIGHT / 2
        self.texture_sprite.width = SCREEN_WIDTH
        self.texture_sprite.height = SCREEN_HEIGHT

        # 2. En Arcade 3.0+ los sprites se dibujan a través de una SpriteList
        self.sprite_list = arcade.SpriteList()
        self.sprite_list.append(self.texture_sprite)

    def on_update(self, delta_time: float):
        # Avanzar el tiempo para darle movimiento a la textura
        self.time_offset += delta_time * 1.5

        # Crear la malla de coordenadas con np.meshgrid() para X e Y
        resolution = 150
        x_vals = np.linspace(-15, 15, resolution)
        y_vals = np.linspace(-15, 15, resolution)
        x, y = np.meshgrid(x_vals, y_vals)

        # Aplicar la fórmula de trigonometría
        matrix = np.sin(x * 0.2 + self.time_offset) + np.cos(y * 0.2 + self.time_offset)

        # Normalizar los valores estrictamente entre 0 y 255 (tipo uint8)
        matrix_min = matrix.min()
        matrix_max = matrix.max()
        normalized = (matrix - matrix_min) / (matrix_max - matrix_min) * 255.0
        img_array = normalized.astype(np.uint8)

        # Convertir con Pillow, forzar RGBA y actualizar la textura del sprite
        pil_image = Image.fromarray(img_array, mode="L").convert("RGBA")
        dynamic_texture = arcade.Texture(image=pil_image)
        self.texture_sprite.texture = dynamic_texture

    def on_draw(self):
        # En Arcade 3.0+ se utiliza self.clear() en lugar de arcade.start_render()
        self.clear()

        # Dibujar la lista de sprites con la textura procedural actualizada
        self.sprite_list.draw()

def main():
    window = ProceduralTextureWindow()
    arcade.run()

if __name__ == "__main__":
    main()
