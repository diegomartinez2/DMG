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

        self.time_offset = 0.0

        # Usamos un sprite para renderizar la textura de forma robusta y moderna
        self.texture_sprite = arcade.Sprite()
        self.texture_sprite.center_x = SCREEN_WIDTH / 2
        self.texture_sprite.center_y = SCREEN_HEIGHT / 2
        self.texture_sprite.width = SCREEN_WIDTH
        self.texture_sprite.height = SCREEN_HEIGHT

    def on_update(self, delta_time: float):
        # Avanzar el tiempo para darle movimiento a la textura
        self.time_offset += delta_time * 1.5

        # 1. Crear la malla de coordenadas con np.meshgrid() para X e Y
        resolution = 150
        x_vals = np.linspace(-15, 15, resolution)
        y_vals = np.linspace(-15, 15, resolution)
        x, y = np.meshgrid(x_vals, y_vals)

        # 2. Aplicar la fórmula de trigonometría
        matrix = np.sin(x * 0.2 + self.time_offset) + np.cos(y * 0.2 + self.time_offset)

        # 3. Normalizar los valores estrictamente entre 0 y 255 (tipo uint8)
        matrix_min = matrix.min()
        matrix_max = matrix.max()
        normalized = (matrix - matrix_min) / (matrix_max - matrix_min) * 255.0
        img_array = normalized.astype(np.uint8)

        # 4. Convertir con Pillow, forzar RGBA y actualizar la textura del sprite
        pil_image = Image.fromarray(img_array, mode="L").convert("RGBA")
        dynamic_texture = arcade.Texture(image=pil_image)
        self.texture_sprite.texture = dynamic_texture

    def on_draw(self):
        arcade.start_render()

        # Dibujar el sprite que contiene nuestra textura procedural dinámica
        self.texture_sprite.draw()

def main():
    window = ProceduralTextureWindow()
    arcade.run()

if __name__ == "__main__":
    main()
