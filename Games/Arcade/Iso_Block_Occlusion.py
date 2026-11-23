#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Iso_Block_Occlusion.py
#  Bloque isométrico elevado, colisión en rejilla y oclusión con Q, A, O, P en Arcade v3.x.
#

import arcade
from PIL import Image, ImageDraw

SCREEN_WIDTH = 800
SCREEN_HEIGHT = 600
SCREEN_TITLE = "Isometría - Bloque Elevado y Oclusión (Q, A, O, P)"

GRID_SIZE = 8
TILE_W = 64
TILE_H = 32
BLOCK_HEIGHT = 36

BLOCK_GX = 4
BLOCK_GY = 4


def make_tile_texture(width: int, height: int, fill_rgba, border_rgba) -> arcade.Texture:
    """Genera en memoria la textura de un rombo isométrico."""
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


def make_elevated_block_texture(width: int, height: int, z_h: int) -> arcade.Texture:
    """
    Genera en memoria un bloque isométrico elevado:
    - Caras laterales visibles abajo (izquierda y derecha).
    - Cara superior en rombo desplazada 'z_h' píxeles hacia arriba en el lienzo PIL.
    """
    total_h = height + z_h
    img = Image.new("RGBA", (width, total_h), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    # En PIL, y=0 es arriba del lienzo; la base del rombo va de y=z_h a y=z_h+height
    base_y = z_h
    half_w = width // 2
    half_h = height // 2

    top_pt = (half_w, base_y)
    right_pt = (width - 1, base_y + half_h)
    bottom_pt = (half_w, base_y + height - 1)
    left_pt = (0, base_y + half_h)

    # Puntos inferiores proyectados hacia abajo en PIL (que es +z en pantalla)
    left_bot = (left_pt[0], left_pt[1] + z_h)
    bottom_bot = (bottom_pt[0], bottom_pt[1] + z_h)
    right_bot = (right_pt[0], right_pt[1] + z_h)

    # Cara lateral izquierda (más oscura)
    left_face = [left_pt, bottom_pt, bottom_bot, left_bot]
    draw.polygon(left_face, fill=(75, 92, 112, 255), outline=(30, 38, 48, 255))

    # Cara lateral derecha (tono intermedio)
    right_face = [bottom_pt, right_pt, right_bot, bottom_bot]
    draw.polygon(right_face, fill=(98, 120, 145, 255), outline=(30, 38, 48, 255))

    # Cara superior (más clara)
    top_face = [top_pt, right_pt, bottom_pt, left_pt]
    draw.polygon(top_face, fill=(135, 165, 198, 255), outline=(220, 232, 245, 255))

    return arcade.Texture(image=img)


def make_sphere_texture(diameter: int, fill_rgba, rim_rgba) -> arcade.Texture:
    """Genera en memoria una esfera sencilla para usar en SpriteList."""
    img = Image.new("RGBA", (diameter, diameter), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    draw.ellipse([1, 1, diameter - 2, diameter - 2], fill=fill_rgba, outline=rim_rgba, width=2)
    return arcade.Texture(image=img)


class IsometricOcclusionWindow(arcade.Window):
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        self.background_color = arcade.color.DARK_SLATE_GRAY

        self.origin_x = SCREEN_WIDTH // 2
        self.origin_y = 130

        # Posición inicial de la esfera (justo detrás del bloque para probar oclusión al arrancar)
        self.gx = 5
        self.gy = 5

        self.tile_texture = make_tile_texture(
            TILE_W, TILE_H,
            fill_rgba=(45, 58, 70, 220),
            border_rgba=(160, 175, 190, 255),
        )
        self.block_texture = make_elevated_block_texture(TILE_W, TILE_H, BLOCK_HEIGHT)
        self.sphere_texture = make_sphere_texture(
            26,
            fill_rgba=(220, 55, 65, 255),
            rim_rgba=(255, 235, 235, 255),
        )

        # Suelo: lo dibujamos primero completo para que sirva de base limpia
        self.floor_list = arcade.SpriteList()
        for gx in range(GRID_SIZE):
            for gy in range(GRID_SIZE):
                sx, sy = self.grid_to_screen(gx, gy)
                sprite = arcade.Sprite(self.tile_texture)
                sprite.center_x = sx
                sprite.center_y = sy
                self.floor_list.append(sprite)

        # Bloque elevado: su ancla inferior en Arcade es la base del rombo (sy),
        # y su textura crece hacia arriba en +BLOCK_HEIGHT/2
        bx, by = self.grid_to_screen(BLOCK_GX, BLOCK_GY)
        self.block_sprite = arcade.Sprite(self.block_texture)
        self.block_sprite.center_x = bx
        self.block_sprite.center_y = by + (BLOCK_HEIGHT / 2.0)

        # Esfera
        self.sphere_sprite = arcade.Sprite(self.sphere_texture)
        self.sync_sphere_position()

        # Lista combinada para ordenar por profundidad (de atrás hacia adelante)
        self.render_order = arcade.SpriteList()
        self.rebuild_render_order()

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

    def rebuild_render_order(self):
        """
        Ordena el bloque y la esfera de atrás hacia adelante:
        Cuanto mayor es (gx + gy), más atrás (arriba) está en pantalla.
        Si empatan en profundidad, dibujamos primero la esfera y luego el bloque
        para que el frente del bloque tape la parte trasera de la esfera.
        """
        self.render_order.clear()
        items = [
            (BLOCK_GX + BLOCK_GY, 1, self.block_sprite),
            (self.gx + self.gy, 0, self.sphere_sprite),
        ]
        # Orden descendente por profundidad (mayor depth primero = más lejos/arriba)
        items.sort(key=lambda entry: (entry[0], entry[1]), reverse=True)
        for _, _, spr in items:
            self.render_order.append(spr)

    def on_key_press(self, symbol: int, modifiers: int):
        """Un paso discreto en la rejilla por cada pulsación de Q, A, O, P."""
        new_gx, new_gy = self.gx, self.gy

        if symbol == arcade.key.Q:
            new_gx += 1
        elif symbol == arcade.key.A:
            new_gx -= 1
        elif symbol == arcade.key.O:
            new_gy += 1
        elif symbol == arcade.key.P:
            new_gy -= 1
        else:
            return

        # Límites del mapa y colisión sólida con el bloque elevado
        if not (0 <= new_gx < GRID_SIZE and 0 <= new_gy < GRID_SIZE):
            return
        if new_gx == BLOCK_GX and new_gy == BLOCK_GY:
            return

        self.gx = new_gx
        self.gy = new_gy
        self.sync_sphere_position()
        self.rebuild_render_order()

    def on_draw(self):
        self.clear()
        self.floor_list.draw()
        self.render_order.draw()

        arcade.draw_text(
            f"Controles: [Q/A] [O/P]  |  Esfera: ({self.gx}, {self.gy})  |  Bloque: ({BLOCK_GX}, {BLOCK_GY})",
            16, SCREEN_HEIGHT - 32, arcade.color.LIGHT_GRAY, 14
        )


def main():
    window = IsometricOcclusionWindow()
    arcade.run()


if __name__ == "__main__":
    main()
