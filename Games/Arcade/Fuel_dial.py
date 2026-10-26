#!/usr/bin/env python
'''
Usa la librería Arcade para crear un indicador analógico.
Actualizado para compatibilidad total con Arcade v3.x (uso de SpriteList).
'''
import math
import arcade

# Configuración de la ventana
SCREEN_WIDTH = 600
SCREEN_HEIGHT = 600
SCREEN_TITLE = "Indicador de Combustible - Arcade"

class FuelGaugeWindow(arcade.Window):
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        self.background_color = arcade.color.DARK_SLATE_GRAY

        # Coordenadas centrales del dial
        self.dial_center_x = SCREEN_WIDTH // 2
        self.dial_center_y = SCREEN_HEIGHT // 2

        # Dimensiones del dial y cálculo del radio
        self.dial_width = 300
        self.dial_height = 300
        self.radius = self.dial_width / 2

        # 1. Crear un SpriteList y el Sprite para el dial (Requerido en Arcade 3.x)
        self.dial_list = arcade.SpriteList()
        try:
            self.dial_sprite = arcade.Sprite("fuel.png")
            self.dial_sprite.center_x = self.dial_center_x
            self.dial_sprite.center_y = self.dial_center_y
            self.dial_sprite.width = self.dial_width
            self.dial_sprite.height = self.dial_height
            self.dial_list.append(self.dial_sprite)
        except Exception:
            self.dial_sprite = None

        # Longitud de la aguja: 2/3 del radio de la imagen
        self.needle_length = self.radius * (2 / 3)

        # Control del tiempo para la animación (en segundos)
        self.time_elapsed = 0.0
        self.duration = 20.0  # Duración del giro de 180º

    def on_update(self, delta_time: float):
        # Acumular el tiempo transcurrido en cada frame
        self.time_elapsed += delta_time

    def on_draw(self):
        # Limpiar la pantalla
        self.clear()

        # 2. Dibujar el dial de combustible usando el SpriteList
        if self.dial_sprite:
            self.dial_list.draw()
        else:
            # Respaldo visual si falta la imagen PNG 'fuel.png'
            arcade.draw_circle_filled(self.dial_center_x, self.dial_center_y, self.radius, arcade.color.LIGHT_GRAY)
            arcade.draw_circle_outline(self.dial_center_x, self.dial_center_y, self.radius, arcade.color.BLACK, 3)

        # 3. Calcular la rotación en sentido antihorario (de 0 a 1 en 20 segundos)
        progress = (self.time_elapsed % self.duration) / self.duration

        # 90° (arriba) a 270° (abajo) en sentido antihorario sumando grados
        angle_deg = 90.0 + (progress * 180.0)
        angle_rad = math.radians(angle_deg)

        # 4. Calcular el extremo de la línea de la aguja
        end_x = self.dial_center_x + self.needle_length * math.cos(angle_rad)
        end_y = self.dial_center_y + self.needle_length * math.sin(angle_rad)

        # Dibujar la aguja (línea roja)
        arcade.draw_line(
            self.dial_center_x, self.dial_center_y,
            end_x, end_y,
            arcade.color.RED,
            line_width=4
        )

        # Perno central de la aguja
        arcade.draw_circle_filled(self.dial_center_x, self.dial_center_y, 6, arcade.color.DARK_RED)

def main():
    window = FuelGaugeWindow()
    arcade.run()

if __name__ == "__main__":
    main()
