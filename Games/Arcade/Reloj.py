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
        """Configura los elementos del juego. Se llama una sola vez al inicio."""

        self.sprite_list = arcade.SpriteList()

        # 1. Cargar el fondo del reloj (512x512) centrado
        self.clock_background = arcade.Sprite(
            "reloj_fondo.png", center_x=CLOCK_CENTER_X, center_y=CLOCK_CENTER_Y
        )

        # 2. Cargar la manecilla de las horas (256x32)
        # Cargamos el sprite de forma nativa sin pasarle el centro en el constructor
        self.hour_hand = arcade.Sprite("manecilla_horas.png")

        # 3. Cargar la manecilla de los minutos (256x32)
        self.minute_hand = arcade.Sprite("manecilla_minutos.png")

        # Añadimos los sprites a la lista en orden de capas
        self.sprite_list.append(self.clock_background)
        self.sprite_list.append(self.hour_hand)
        self.sprite_list.append(self.minute_hand)

        # --- INTERFAZ GRÁFICA (GUI) ---
        v_box = arcade.gui.UIBoxLayout()
        exit_button = arcade.gui.UIFlatButton(text="Salir del Juego", width=150)
        v_box.add(exit_button)

        @exit_button.event("on_click")
        def on_click_exit(event):
            arcade.exit()

        anchor_layout = arcade.gui.UIAnchorLayout()
        anchor_layout.add(
            child=v_box,
            anchor_x="left",
            anchor_y="bottom",
            align_x=20,
            align_y=20,
        )
        self.gui_manager.add(anchor_layout)

    def on_update(self, delta_time):
        """Lógica del juego: Se ejecuta automáticamente unas 60 veces por segundo."""

        # Obtenemos la hora actual del sistema operativo
        current_time = time.localtime()
        hours = current_time.tm_hour % 12
        minutes = current_time.tm_min
        seconds = current_time.tm_sec

        # CÁLCULO DE ÁNGULOS EN SENTIDO HORARIO:
        # Partimos de 90° (las 12 en punto) y restamos para avanzar hacia la derecha.
        minute_angle = 90 - (minutes * 6) - (seconds * 0.1)
        hour_angle = 90 - (hours * 30) - (minutes * 0.5)

        # Aplicamos la rotación a los objetos
        self.minute_hand.angle = minute_angle
        self.hour_hand.angle = hour_angle

        # SOLUCIÓN DE ANCLAJE COMPATIBLE CON ARCADE 3.x:
        # Como el motor rota sobre el centro del sprite de forma fija, recalculamos la posición
        # del centro del sprite basándonos en el ángulo actual para mantener el píxel 16 estático.
        # La distancia desde el centro del sprite (128) al píxel de anclaje (16) es de 112 píxeles.
        distancia_eje = 128 - 16

        # Convertimos los ángulos a radianes para la librería math
        rad_min = math.radians(minute_angle)
        rad_hour = math.radians(hour_angle)

        # Reposicionamos los centros de los sprites para contrarrestar la rotación
        self.minute_hand.center_x = CLOCK_CENTER_X + distancia_eje * math.cos(rad_min)
        self.minute_hand.center_y = CLOCK_CENTER_Y + distancia_eje * math.sin(rad_min)

        self.hour_hand.center_x = CLOCK_CENTER_X + distancia_eje * math.cos(rad_hour)
        self.hour_hand.center_y = CLOCK_CENTER_Y + distancia_eje * math.sin(rad_hour)

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
