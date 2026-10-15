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

        # 1. Cargar la imagen del dial (asegúrate de tener un archivo 'fuel_dial.png' con transparencia en la misma carpeta)
        try:
            self.dial_texture = arcade.load_texture("fuel.png")
        except Exception:
            # Si no se encuentra la imagen, se manejará de forma gráfica alternativa en el dibujado
            self.dial_texture = None

        # Coordenadas del centro de la pantalla (donde se ubicará el dial)
        self.center_x = SCREEN_WIDTH // 2
        self.center_y = SCREEN_HEIGHT // 2

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

        # 2. Dibujar el dial de combustible (imagen PNG con transparencia)
        if self.dial_texture:
            arcade.draw_texture_rectangle(
                self.center_x, self.center_y,
                self.dial_width, self.dial_height,
                self.dial_texture
            )
        else:
            # Respaldo visual por si falta la imagen PNG en local
            arcade.draw_circle_filled(self.center_x, self.center_y, self.radius, arcade.color.LIGHT_GRAY)
            arcade.draw_circle_outline(self.center_x, self.center_y, self.radius, arcade.color.BLACK, 3)

        # 3. Calcular la rotación en sentido antihorario (de 0 a 1 en 20 segundos, de forma cíclica)
        progress = (self.time_elapsed % self.duration) / self.duration

        # En trigonometría estándar (y en Arcade/OpenGL):
        # 90 grados = Arriba
        # 180 grados = Izquierda
        # 270 grados = Abajo
        # Moverse de 90° a 270° incrementando el ángulo realiza un giro antihorario (180º exactos).
        angle_deg = 90.0 + (progress * 180.0)
        angle_rad = math.radians(angle_deg)

        # 4. Calcular el extremo de la línea roja usando trigonometría
        end_x = self.center_x + self.needle_length * math.cos(angle_rad)
        end_y = self.center_y + self.needle_length * math.sin(angle_rad)

        # Dibujar la aguja (línea roja) desde el centro hacia el extremo calculado
        arcade.draw_line(
            self.center_x, self.center_y,
            end_x, end_y,
            arcade.color.RED,
            line_width=4
        )

        # Opcional: un pequeño perno central oscuro para mejorar el aspecto visual de la aguja
        arcade.draw_circle_filled(self.center_x, self.center_y, 6, arcade.color.DARK_RED)

def main():
    window = FuelGaugeWindow()
    arcade.run()

if __name__ == "__main__":
    main()
