#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Iso_Jump_Sine.py
#  Salto isométrico con curva seno (misma altura, subida y bajada) en Arcade v3.x.
#

import math
import arcade
from PIL import Image, ImageDraw

SCREEN_WIDTH = 800
SCREEN_HEIGHT = 600
SCREEN_TITLE = "Isometría - Salto Senoidal y Alturas (Q, A, O, P + ESPACIO)"

GRID_SIZE = 8
TILE_W = 64
TILE_H = 32
BLOCK_HEIGHT = 40

BLOCK_GX = 4
BLOCK_GY = 4

JUMP_DURATION = 0.42     # Duración del salto en segundos
JUMP_PEAK_EXTRA = 28.0   # Altura extra por encima del punto más alto entre inicio y fin
MAX_STEP_UP = BLOCK_HEIGHT # Máxima diferencia de altura que puede subir un salto


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
    """Genera en memoria un bloque isométrico elevado con altura z_h."""
    total_h = height + z_h
    img = Image.new("RGBA", (width, total_h), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    half_w = width // 2
    half_h = height // 2

    top_pt = (half_w, 0)
    right_pt = (width - 1, half_h)
    bottom_pt = (half_w, height - 1)
    left_pt = (0, half_h)

    left_bot = (left_pt[0], min(total_h - 1, left_pt[1] + z_h))
    bottom_bot = (bottom_pt[0], min(total_h - 1, bottom_pt[1] + z_h))
    right_bot = (right_pt[0], min(total_h - 1, right_pt[1] + z_h))

    draw.polygon([left_pt, bottom_pt, bottom_bot, left_bot], fill=(75, 92, 112, 255), outline=(30, 38, 48, 255))
    draw.polygon([bottom_pt, right_pt, right_bot, bottom_bot], fill=(98, 120, 145, 255), outline=(30, 38, 48, 255))
    draw.polygon([top_pt, right_pt, bottom_pt, left_pt], fill=(135, 165, 198, 255), outline=(220, 232, 245, 255))

    return arcade.Texture(image=img)


def make_sphere_texture(diameter: int, fill_rgba, rim_rgba) -> arcade.Texture:
    """Genera en memoria una esfera sencilla para usar en SpriteList."""
    img = Image.new("RGBA", (diameter, diameter), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    draw.ellipse([1, 1, diameter - 2, diameter - 2], fill=fill_rgba, outline=rim_rgba, width=2)
    return arcade.Texture(image=img)


class IsometricJumpWindow(arcade.Window):
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        self.background_color = arcade.color.DARK_SLATE_GRAY

        self.origin_x = SCREEN_WIDTH // 2
        self.origin_y = 130

        # Posición lógica y altura actual de apoyo de la esfera
        self.gx = 5.0
        self.gy = 4.0
        self.z_current = 0.0
        self.last_dir = (-1, 0)  # Hacia (4, 4) donde está el bloque

        # Estado de salto (None si está en reposo)
        self.jump_data = None
        self.status_msg = "Listo: usa Q/A/O/P para moverte y ESPACIO para saltar."

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

        # Suelo base
        self.floor_list = arcade.SpriteList()
        for gx in range(GRID_SIZE):
            for gy in range(GRID_SIZE):
                sx, sy = self.grid_to_screen(gx, gy)
                sprite = arcade.Sprite(self.tile_texture)
                sprite.center_x = sx
                sprite.center_y = sy
                self.floor_list.append(sprite)

        # Bloque elevado
        bx, by = self.grid_to_screen(BLOCK_GX, BLOCK_GY)
        self.block_sprite = arcade.Sprite(self.block_texture)
        self.block_sprite.center_x = bx
        self.block_sprite.center_y = by + (BLOCK_HEIGHT / 2.0)

        # Esfera
        self.sphere_sprite = arcade.Sprite(self.sphere_texture)
        self.render_order = arcade.SpriteList()

        self.sync_sphere_visuals()
        self.rebuild_render_order()

    def cell_elevation(self, gx: int, gy: int) -> float:
        """Devuelve la altura Z del suelo/superficie en la casilla (gx, gy)."""
        if int(round(gx)) == BLOCK_GX and int(round(gy)) == BLOCK_GY:
            return float(BLOCK_HEIGHT)
        return 0.0

    def grid_to_screen(self, gx: float, gy: float, z_height: float = 0.0):
        """Convierte coordenadas lógicas (gx, gy) y altura z a coordenadas de pantalla."""
        sx = self.origin_x + (gx - gy) * (TILE_W // 2)
        sy = self.origin_y + (gx + gy) * (TILE_H // 2) + z_height
        return sx, sy

    def sync_sphere_visuals(self):
        """Actualiza la posición en pantalla de la esfera según (gx, gy) y z_current."""
        sx, sy = self.grid_to_screen(self.gx, self.gy, z_height=self.z_current + 10.0)
        self.sphere_sprite.center_x = sx
        self.sphere_sprite.center_y = sy

    def rebuild_render_order(self):
        """
        Ordena el bloque y la esfera de atrás hacia adelante:
        - Si la esfera está por encima de la cima del bloque, va después (encima del bloque).
        - En caso contrario, ordena por profundidad isométrica (gx + gy descendente).
        """
        self.render_order.clear()
        sphere_depth = self.gx + self.gy
        block_depth = BLOCK_GX + BLOCK_GY

        if self.z_current >= BLOCK_HEIGHT - 0.5:
            # Está sobre el bloque: dibujamos primero el bloque y luego la esfera encima
            self.render_order.append(self.block_sprite)
            self.render_order.append(self.sphere_sprite)
        else:
            items = [
                (block_depth, 1, self.block_sprite),
                (sphere_depth, 0, self.sphere_sprite),
            ]
            items.sort(key=lambda entry: (entry[0], entry[1]), reverse=True)
            for _, _, spr in items:
                self.render_order.append(spr)

    def try_start_jump(self):
        """Inicia un salto de 1 casilla en la última dirección pulsada si es válido."""
        if self.jump_data is not None:
            return

        dx, dy = self.last_dir
        start_gx, start_gy = self.gx, self.gy
        end_gx, end_gy = start_gx + dx, start_gy + dy

        # Comprobar límites de la rejilla
        if not (0 <= end_gx < GRID_SIZE and 0 <= end_gy < GRID_SIZE):
            self.status_msg = "Fuera del tablero: no puedes saltar ahí."
            return

        z_start = self.cell_elevation(start_gx, start_gy)
        z_end = self.cell_elevation(end_gx, end_gy)

        # Comprobar límite de subida
        if z_end - z_start > MAX_STEP_UP:
            self.status_msg = "Demasiado alto para saltar de una vez."
            return

        # Amplitud del arco seno: supera el desnivel positivo para que el salto sea nítido
        delta_z = z_end - z_start
        amp = JUMP_PEAK_EXTRA + max(0.0, delta_z)

        if delta_z > 0:
            tipo = "SUBIDA"
        elif delta_z < 0:
            tipo = "BAJADA"
        else:
            tipo = "MISMA ALTURA"

        self.status_msg = f"Saltando ({tipo}): Z {int(z_start)} -> {int(z_end)}"
        self.jump_data = {
            "elapsed": 0.0,
            "duration": JUMP_DURATION,
            "start_gx": start_gx,
            "start_gy": start_gy,
            "end_gx": end_gx,
            "end_gy": end_gy,
            "z_start": z_start,
            "z_end": z_end,
            "amp": amp,
        }

    def on_key_press(self, symbol: int, modifiers: int):
        """Movimiento paso a paso con Q/A/O/P y salto con ESPACIO."""
        if symbol == arcade.key.SPACE:
            self.try_start_jump()
            return

        if self.jump_data is not None:
            return

        dx, dy = 0, 0
        if symbol == arcade.key.Q:
            dx = 1
        elif symbol == arcade.key.A:
            dx = -1
        elif symbol == arcade.key.O:
            dy = 1
        elif symbol == arcade.key.P:
            dy = -1
        else:
            return

        self.last_dir = (dx, dy)
        new_gx = self.gx + dx
        new_gy = self.gy + dy

        if not (0 <= new_gx < GRID_SIZE and 0 <= new_gy < GRID_SIZE):
            return

        target_z = self.cell_elevation(new_gx, new_gy)
        # Un paso normal no sube a un bloque elevado; usa ESPACIO para saltar
        if target_z > self.z_current:
            self.status_msg = f"Bloque elevado en ({new_gx}, {new_gy}): pulsa ESPACIO para saltar."
            return

        self.gx = new_gx
        self.gy = new_gy
        self.z_current = target_z
        self.status_msg = f"Posición ({int(self.gx)}, {int(self.gy)}) | Z={int(self.z_current)}"
        self.sync_sphere_visuals()
        self.rebuild_render_order()

    def on_update(self, delta_time: float):
        if self.jump_data is None:
            return

        jd = self.jump_data
        jd["elapsed"] += min(delta_time, 0.05)
        u = min(1.0, jd["elapsed"] / jd["duration"])

        # Interpolación horizontal en la rejilla
        self.gx = jd["start_gx"] + (jd["end_gx"] - jd["start_gx"]) * u
        self.gy = jd["start_gy"] + (jd["end_gy"] - jd["start_gy"]) * u

        # Cota base interpolada + arco sinusoidal
        base_z = jd["z_start"] + (jd["z_end"] - jd["z_start"]) * u
        arc_z = jd["amp"] * math.sin(math.pi * u)
        self.z_current = base_z + arc_z

        if u >= 1.0:
            self.gx = float(jd["end_gx"])
            self.gy = float(jd["end_gy"])
            self.z_current = jd["z_end"]
            self.jump_data = None
            self.status_msg = f"Aterrizaje en ({int(self.gx)}, {int(self.gy)}) | Z={int(self.z_current)}"

        self.sync_sphere_visuals()
        self.rebuild_render_order()

    def on_draw(self):
        self.clear()
        self.floor_list.draw()
        self.render_order.draw()

        arcade.draw_text(
            f"[Q/A/O/P] Mover/Orientar   [ESPACIO] Saltar   |   {self.status_msg}",
            16, SCREEN_HEIGHT - 32, arcade.color.LIGHT_GRAY, 13
        )


def main():
    window = IsometricJumpWindow()
    arcade.run()


if __name__ == "__main__":
    main()
