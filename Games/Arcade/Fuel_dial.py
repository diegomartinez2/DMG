#!/usr/bin/env python
'''
Usa la librería Arcade para crear un indicador analógico.
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
        arcade.set_background_color(arcade.color.DARK_SLATE_GRAY)

        # 1. Cargar la imagen del dial (asegúrate de tener 'fuel_dial.png' en la misma carpeta)
        try:
            self.dial_texture = arcade.load_texture("fuel.png")
        except Exception:
            self.dial_texture = None

        # CORREGIDO: Usar nombres personalizados para evitar colisión con arcade.Window
        self.dial_center_x = SCREEN_WIDTH // 2
        self.dial_center_y = SCREEN_HEIGHT // 2

        # Dimensiones del dial y cálculo del radio
        self.dial_width = 300
        self.dial_height = 300
        self.radius = self.dial_width / 2

        # Longitud de la aguja: 2/3 del radio de la imagen
        self.needle_length = self.radius * (2 / 3)

        # Control del tiempo para la animación (en segundos)
        self.time_elapsed = 0.0
        self.duration = 20.0  # Duración del giro de 180º

    def on_update(self, delta_time: float):
        # Acumular el tiempo transcurrido en cada frame
        self.time_elapsed += delta_time

    def on_draw(self):
        # Iniciar el proceso de renderizado
        arcade.start_render()

        # 2. Dibujar el dial de combustible
        if self.dial_texture:
            arcade.draw_texture_rectangle(
                self.dial_center_x, self.dial_center_y,
                self.dial_width, self.dial_height,
                self.dial_texture
            )
        else:
            # Respaldo visual si falta la imagen PNG
            arcade.draw_circle_filled(self.dial_center_x, self.dial_center_y, self.radius, arcade.color.LIGHT_GRAY)
            arcade.draw_circle_outline(self.dial_center_x, self.dial_center_y, self.radius, arcade.color.BLACK, 3)

        # 3. Calcular la rotación en sentido antihorario (de 0 a 1 en 20 segundos)
        progress = (self.time_elapsed % self.duration) / self.duration

        # 90° (arriba) a 270° (abajo) en sentido antihorario sumando grados
        angle_deg = 90.0 + (progress * 180.0)
        angle_rad = math.radians(angle_deg)

        # 4. Calcular el extremo de la línea roja
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
