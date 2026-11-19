#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Iso_Q_A_O_P.py
#  Movimiento discreto por pasos sobre rejilla isométrica en Arcade v3.x con Q, A, O, P.
#

import arcade
from PIL import Image, ImageDraw

SCREEN_WIDTH = 800
SCREEN_HEIGHT = 600
SCREEN_TITLE = "Isometría Clásica - Control Q, A, O, P"

GRID_SIZE = 8
TILE_W = 64
TILE_H = 32


def make_tile_texture(width: int, height: int, fill_rgba, border_rgba) -> arcade.Texture:
    """Genera en memoria la textura de un rombo isométrico sin tocar el disco."""
    img = Image.new("RGBA", (width, height), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    diamond = [
        (width // 2, 0),
        (width - 1, height // 2),
        (width // 2, height - 1),
        (0, height // 2),
    ]
    draw.polygon(diamond, fill=fill_rgba, outline=border_rgba)
    return arcade.Texture(image=img)


def make_sphere_texture(diameter: int, fill_rgba, rim_rgba) -> arcade.Texture:
    """Genera en memoria una esfera sencilla para usar en SpriteList."""
    img = Image.new("RGBA", (diameter, diameter), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    draw.ellipse([1, 1, diameter - 2, diameter - 2], fill=fill_rgba, outline=rim_rgba, width=2)
    return arcade.Texture(image=img)


class IsometricQaopWindow(arcade.Window):
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        self.background_color = arcade.color.DARK_SLATE_GRAY

        # Origen de la esquina (gx=0, gy=0) en pantalla
        self.origin_x = SCREEN_WIDTH // 2
        self.origin_y = 140

        # Posición lógica discreta de la esfera en la rejilla
        self.gx = 3
        self.gy = 3

        # Texturas procedurales ligeras en memoria
        self.tile_texture = make_tile_texture(
            TILE_W, TILE_H,
            fill_rgba=(45, 58, 70, 220),
            border_rgba=(160, 175, 190, 255),
        )
        self.sphere_texture = make_sphere_texture(
            26,
            fill_rgba=(220, 55, 65, 255),
            rim_rgba=(255, 235, 235, 255),
        )

        # Creamos los sprites de las teselas y los ordenamos de atrás hacia adelante (gx + gy)
        self.tiles_by_depth = []
        for gx in range(GRID_SIZE):
            for gy in range(GRID_SIZE):
                sx, sy = self.grid_to_screen(gx, gy)
                sprite = arcade.Sprite(self.tile_texture)
                sprite.center_x = sx
                sprite.center_y = sy
                depth = gx + gy
                self.tiles_by_depth.append((depth, gx, gy, sprite))

        self.tiles_by_depth.sort(key=lambda item: item[0])

        # Sprite de la esfera
        self.sphere_sprite = arcade.Sprite(self.sphere_texture)
        self.sync_sphere_position()

    def grid_to_screen(self, gx: int, gy: int, z_height: float = 0.0):
        """Convierte coordenadas lógicas (gx, gy) y altura z a coordenadas de pantalla."""
        sx = self.origin_x + (gx - gy) * (TILE_W // 2)
        sy = self.origin_y + (gx + gy) * (TILE_H // 2) + z_height
        return sx, sy

    def sync_sphere_position(self):
        """Coloca la esfera ligeramente elevada sobre el centro de su tesela."""
        sx, sy = self.grid_to_screen(self.gx, self.gy, z_height=10.0)
        self.sphere_sprite.center_x = sx
        self.sphere_sprite.center_y = sy

    def on_key_press(self, symbol: int, modifiers: int):
        """Un paso discreto en la rejilla por cada pulsación de Q, A, O, P."""
        new_gx, new_gy = self.gx, self.gy

        if symbol == arcade.key.Q:
            new_gx -= 1
        elif symbol == arcade.key.A:
            new_gx += 1
        elif symbol == arcade.key.O:
            new_gy += 1
        elif symbol == arcade.key.P:
            new_gy -= 1
        else:
            return

        # Mantener dentro de los límites de la rejilla
        if 0 <= new_gx < GRID_SIZE and 0 <= new_gy < GRID_SIZE:
            self.gx = new_gx
            self.gy = new_gy
            self.sync_sphere_position()

    def on_draw(self):
        self.clear()

        sphere_depth = self.gx + self.gy
        sphere_drawn = False

        # Pintamos primero las teselas con menor profundidad (más atrás);
        # cuando la profundidad alcanza la de la esfera, pintamos la esfera encima
        # y luego las teselas que quedan por delante de ella.
        drawn_sphere_list = arcade.SpriteList()
        drawn_sphere_list.append(self.sphere_sprite)

        front_tiles = arcade.SpriteList()
        back_tiles = arcade.SpriteList()

        for depth, _, _, sprite in self.tiles_by_depth:
            if depth <= sphere_depth:
                back_tiles.append(sprite)
            else:
                front_tiles.append(sprite)

        back_tiles.draw()
        drawn_sphere_list.draw()
        front_tiles.draw()

        arcade.draw_text(
            f"Controles: [Q/A] eje X  [O/P] eje Y   |   Pos: ({self.gx}, {self.gy})",
            16, SCREEN_HEIGHT - 32, arcade.color.LIGHT_GRAY, 14
        )


def main():
    window = IsometricQaopWindow()
    arcade.run()


if __name__ == "__main__":
    main()
