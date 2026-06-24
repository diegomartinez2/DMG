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

        # 1. Cargar el fondo del reloj (512x512)
        self.clock_background = arcade.Sprite(
            "reloj_fondo.png", center_x=CLOCK_CENTER_X, center_y=CLOCK_CENTER_Y
        )

        # 2. Cargar la manecilla de las horas (256x32)
        self.hour_hand = arcade.Sprite("manecilla_horas.png")

        # CORRECCIÓN DE ANCLAJE EN ARCADE 3.x:
        # La imagen mide 256 de ancho. Su centro es 128. Queremos que rote en el píxel 16.
        # Desplazamos el punto de dibujo (center) respecto al eje del reloj para simular el pivote.
        self.hour_hand.center_x = CLOCK_CENTER_X + (128 - 16)
        self.hour_hand.center_y = CLOCK_CENTER_Y

        # 3. Cargar la manecilla de los minutos (256x32) con la misma lógica
        self.minute_hand = arcade.Sprite("manecilla_minutos.png")
        self.minute_hand.center_x = CLOCK_CENTER_X + (128 - 16)
        self.minute_hand.center_y = CLOCK_CENTER_Y

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
        """Lógica del juego: Se ejecuta automáticamente unas 60 veces por second."""

        # Obtenemos la hora actual del sistema operativo
        current_time = time.localtime()
        hours = current_time.tm_hour % 12
        minutes = current_time.tm_min
        seconds = current_time.tm_sec

        # CORRECCIÓN DE SENTIDO HORARIO:
        # En Arcade 0° es a la derecha (las 3 en un reloj).
        # Para que empiece en las 12 (arriba), partimos de 90°.
        # Al usar el signo MENOS (-), forzamos a que el ángulo avance hacia la derecha (sentido horario).

        minute_angle = 90 - (minutes * 6) - (seconds * 0.1)
        hour_angle = 90 - (hours * 30) - (minutes * 0.5)

        # Aplicamos la rotación
        self.minute_hand.angle = minute_angle
        self.hour_hand.angle = hour_angle

        # CORRECCIÓN DINÁMICA DE ROTACIÓN SOBRE PIVOTE EXCENTRICO:
        # Como Arcade 3.x rota los sprites sobre su centro geométrico, si cambiamos el ángulo,
        # debemos recalcular matemáticamente dónde queda el centro del sprite para que el extremo (el píxel 16)
        # se quede "clavado" inmóvil en el centro del reloj (CLOCK_CENTER_X, CLOCK_CENTER_Y).

        distancia_al_pivote = 128 - 16  # Distancia desde el centro del sprite al píxel del eje

        # Convertimos el ángulo a radianes para las funciones de math
        rad_min = math.radians(minute_angle)
        rad_hour = math.radians(hour_angle)

        # Aplicamos trigonometría (Órbita circular del centro del sprite alrededor del eje del reloj)
        self.minute_hand.center_x = CLOCK_CENTER_X + distancia_al_pivote * math.cos(rad_min)
        self.minute_hand.center_y = CLOCK_CENTER_Y + distancia_al_pivote * math.sin(rad_min)

        self.hour_hand.center_x = CLOCK_CENTER_X + distancia_al_pivote * math.cos(rad_hour)
        self.hour_hand.center_y = CLOCK_CENTER_Y + distancia_al_pivote * math.sin(rad_hour)

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
