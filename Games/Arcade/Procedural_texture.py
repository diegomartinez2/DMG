import numpy as np
import arcade
from PIL import Image

SCREEN_WIDTH = 600
SCREEN_HEIGHT = 600
SCREEN_TITLE = "Textura Dinámica con NumPy y Arcade"

class ProceduralTextureWindow(arcade.Window):
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE)
        arcade.set_background_color(arcade.color.BLACK)

        # Variable de tiempo para animar el patrón trigonométrico
        self.time_offset = 0.0
        self.dynamic_texture = None

    def on_update(self, delta_time: float):
        # Avanzar el tiempo para darle movimiento a la textura
        self.time_offset += delta_time * 1.5

        # 1. Crear la malla de coordenadas con np.meshgrid() para X e Y
        resolution = 150  # Resolución de la matriz (píxeles lógicos del patrón)
        x_vals = np.linspace(-15, 15, resolution)
        y_vals = np.linspace(-15, 15, resolution)
        x, y = np.meshgrid(x_vals, y_vals)

        # 2. Aplicar la fórmula de trigonometría (añadimos time_offset para dinamismo)
        matrix = np.sin(x * 0.2 + self.time_offset) + np.cos(y * 0.2 + self.time_offset)

        # 3. Normalizar los valores estrictamente entre 0 y 255 (tipo uint8 para escala de grises)
        matrix_min = matrix.min()
        matrix_max = matrix.max()
        normalized = (matrix - matrix_min) / (matrix_max - matrix_min) * 255.0
        img_array = normalized.astype(np.uint8)

        # 4. Convertir la matriz en imagen con Pillow y envolverla como textura en Arcade
        pil_image = Image.fromarray(img_array, mode="L")  # "L" indica escala de grises (luminancia)

        # Inyectamos la imagen como textura dinámica en Arcade con un nombre único
        self.dynamic_texture = arcade.Texture(name="procedural_pattern", image=pil_image)

    def on_draw(self):
        arcade.start_render()

        # Dibujar la textura generada proceduralmente cubriendo toda la ventana
        if self.dynamic_texture:
            arcade.draw_texture_rectangle(
                center_x=SCREEN_WIDTH / 2,
                center_y=SCREEN_HEIGHT / 2,
                width=SCREEN_WIDTH,
                height=SCREEN_HEIGHT,
                texture=self.dynamic_texture
            )

def main():
    window = ProceduralTextureWindow()
    arcade.run()

if __name__ == "__main__":
    main()
