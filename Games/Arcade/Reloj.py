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
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#   MA 02110-1301, USA.
#
#
"""
Ejemplo del uso de las librerias arcade para hacer un juego.
"""
import math
import time
import arcade
import arcade.gui

# --- CONSTANTES ---
SCREEN_WIDTH = 600
SCREEN_HEIGHT = 600
SCREEN_TITLE = "Reloj Analógico con Python Arcade"

# Coordenadas del centro de la pantalla donde irá el reloj
CLOCK_CENTER_X = SCREEN_WIDTH // 2
CLOCK_CENTER_Y = SCREEN_HEIGHT // 2


class ClockGame(arcade.Window):
    """Clase principal del juego que hereda de arcade.Window"""

    def __init__(self):
        # Inicializa la ventana del juego con el tamaño y título
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)

        # Establece el fondo de la ventana (gris oscuro)
        arcade.set_background_color(arcade.color.DARK_SLATE_GRAY)

        # SpriteList para renderizar objetos en Arcade 3.x
        self.sprite_list = None

        # Declaración de variables para los Sprites (imágenes)
        self.clock_background = None
        self.hour_hand = None
        self.minute_hand = None

        # Administrador para la Interfaz Gráfica de Usuario (GUI)
        self.gui_manager = arcade.gui.UIManager()
        self.gui_manager.enable()

    def setup(self):
        """Configura los elementos del juego."""
        self.sprite_list = arcade.SpriteList()

        # Al estar el eje en el centro de las imágenes, solo las colocas en el centro del reloj
        self.clock_background = arcade.Sprite("reloj_fondo.png", center_x=CLOCK_CENTER_X, center_y=CLOCK_CENTER_Y)

        # Las nuevas manecillas extendidas con transparencia
        self.hour_hand = arcade.Sprite("manecilla_horas_nueva.png", center_x=CLOCK_CENTER_X, center_y=CLOCK_CENTER_Y)
        self.minute_hand = arcade.Sprite("manecilla_minutos_nueva.png", center_x=CLOCK_CENTER_X, center_y=CLOCK_CENTER_Y)

        self.sprite_list.append(self.clock_background)
        self.sprite_list.append(self.hour_hand)
        self.sprite_list.append(self.minute_hand)

    def on_update(self, delta_time):
        """Lógica del juego."""
        current_time = time.localtime()
        hours = current_time.tm_hour % 12
        minutes = current_time.tm_min
        seconds = current_time.tm_sec

        # Calculas los ángulos en sentido horario
        minute_angle = 90 - (minutes * 6) - (seconds * 0.1)
        hour_angle = 90 - (hours * 30) - (minutes * 0.5)

        # Al cambiar el ángulo, rotarán perfectamente sobre su propio centro (el eje real)
        self.minute_hand.angle = minute_angle
        self.hour_hand.angle = hour_angle

    def on_draw(self):
        """Renderizado: Dibuja todo en la pantalla."""
        self.clear()
        self.sprite_list.draw()
        self.gui_manager.draw()


def main():
    window = ClockGame()
    window.setup()
    arcade.run()


if __name__ == "__main__":
    main()
