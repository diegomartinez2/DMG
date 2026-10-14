#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Reloj.py
#
#   Copyright 2026 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
"""
Ejemplo de reloj analógico optimizado usando escala de Sprites en Arcade 3.x.
Con esto se puede aprender los canceptos básicos de la librería Arcade.
"""
import time
import arcade
import arcade.gui

# --- CONSTANTES ---
SCREEN_WIDTH = 600
SCREEN_HEIGHT = 600
SCREEN_TITLE = "Reloj Optimizado - Aprendiendo Escalado"

CLOCK_CENTER_X = SCREEN_WIDTH // 2
CLOCK_CENTER_Y = SCREEN_HEIGHT // 2


class ClockGame(arcade.Window):
    """Clase principal del juego"""

    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        arcade.set_background_color(arcade.color.DARK_SLATE_GRAY)

        # Contenedor optimizado para los Sprites
        self.sprite_list = None

        # Variables para los Sprites
        self.clock_background = None
        self.hour_hand = None
        self.minute_hand = None

        # Administrador de la interfaz gráfica (GUI)
        self.gui_manager = arcade.gui.UIManager()
        self.gui_manager.enable()

    def setup(self):
        """Configura los elementos del juego. Se llama una sola vez."""
        self.sprite_list = arcade.SpriteList()

        # 1. Cargar el fondo del reloj
        self.clock_background = arcade.Sprite(
            "reloj_fondo.png", center_x=CLOCK_CENTER_X, center_y=CLOCK_CENTER_Y
        )

        # 2. MANECILLA DE LOS MINUTOS (Tamaño original, escala = 1.0)
        self.minute_hand = arcade.Sprite(
            "manecilla.png",
            center_x=CLOCK_CENTER_X,
            center_y=CLOCK_CENTER_Y,
            scale=1.0
        )

        # 3. MANECILLA DE LAS HORAS (Reutiliza la misma imagen, pero más corta y gruesa)
        # Al usar 'scale', Arcade reduce proporcionalmente el ancho y el alto.
        # Una escala de 0.65 significa que medirá el 65% de la original.
        self.hour_hand = arcade.Sprite(
            "manecilla.png",
            center_x=CLOCK_CENTER_X,
            center_y=CLOCK_CENTER_Y,
            scale=0.65
        )

        # Re-escalamos para hacer la agujas del reloj de tamaños adecuados.
        self.hour_hand.scale_y = 0.80  # s_x
        self.minute_hand.scale_y = 0.50  # s_y

        # Añadimos los elementos en orden de capas
        self.sprite_list.append(self.clock_background)
        self.sprite_list.append(self.hour_hand)       # Las horas abajo
        self.sprite_list.append(self.minute_hand)     # Los minutos encima

        # --- INTERFAZ GRÁFICA (Botón Salir) ---
        v_box = arcade.gui.UIBoxLayout()
        exit_button = arcade.gui.UIFlatButton(text="Salir del Juego", width=150)
        v_box.add(exit_button)

        @exit_button.event("on_click")
        def on_click_exit(event):
            arcade.exit()

        anchor_layout = arcade.gui.UIAnchorLayout()
        anchor_layout.add(
            child=v_box, anchor_x="left", anchor_y="bottom", align_x=20, align_y=20
        )
        self.gui_manager.add(anchor_layout)

    def on_update(self, delta_time):
        """Lógica del juego: Actualiza los ángulos unas 60 veces por segundo."""
        current_time = time.localtime()
        hours = current_time.tm_hour % 12
        minutes = current_time.tm_min
        seconds = current_time.tm_sec

        # Cálculo de ángulos en sentido horario (partiendo de las 12 en punto = 90°)
        minute_angle = 90 + (minutes * 6) + (seconds * 0.1)
        hour_angle = 90 + (hours * 30) + (minutes * 0.5)

        # Aplicamos la rotación. Al estar el eje en el centro perfecto de la imagen,
        # el motor gráfico las hace girar sin descolocarse lo más mínimo.
        self.minute_hand.angle = minute_angle
        self.hour_hand.angle = hour_angle

    def on_draw(self):
        """Renderizado de la pantalla."""
        self.clear()
        self.sprite_list.draw()
        self.gui_manager.draw()


def main():
    window = ClockGame()
    window.setup()
    arcade.run()


if __name__ == "__main__":
    main()
