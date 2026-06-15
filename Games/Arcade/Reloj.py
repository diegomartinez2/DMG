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

        # Declaración de variables para los Sprites (imágenes)
        self.clock_background = None
        self.hour_hand = None
        self.minute_hand = None

        # Variables para almacenar los ángulos de rotación
        self.hour_angle = 0.0
        self.minute_angle = 0.0

        # Administrador para la Interfaz Gráfica de Usuario (GUI)
        self.gui_manager = arcade.gui.UIManager()
        self.gui_manager.enable()

    def setup(self):
        """Configura los elementos del juego. Se llama una sola vez al inicio."""

        # 1. Cargar el fondo del reloj (512x512) en el centro de la pantalla
        self.clock_background = arcade.Sprite(
            "reloj_fondo.png", center_x=CLOCK_CENTER_X, center_y=CLOCK_CENTER_Y
        )

        # 2. Cargar la manecilla de las horas (256x32)
        # Usamos 'anchor_x' para decirle a Arcade que su eje de rotación está a 16 píxeles de la izquierda
        self.hour_hand = arcade.Sprite(
            "manecilla_horas.png",
            center_x=CLOCK_CENTER_X,
            center_y=CLOCK_CENTER_Y,
            anchor_x=16,
            anchor_y=16,
        )

        # 3. Cargar la manecilla de los minutos (256x32) con el mismo punto de anclaje
        self.minute_hand = arcade.Sprite(
            "manecilla_minutos.png",
            center_x=CLOCK_CENTER_X,
            center_y=CLOCK_CENTER_Y,
            anchor_x=16,
            anchor_y=16,
        )

        # --- EJEMPLO DE ELEMENTO GUI (Botón de Salida) ---
        # Creamos una caja vertical para organizar los elementos de la interfaz
        v_box = arcade.gui.UIBoxLayout()

        # Creamos un botón usando arcade.gui
        exit_button = arcade.gui.UIFlatButton(text="Salir del Juego", width=150)
        v_box.add(exit_button)

        # Programamos el evento: qué pasa cuando se hace clic en el botón
        @exit_button.event("on_click")
        def on_click_exit(event):
            arcade.exit()

        # Añadimos la caja de la GUI a la pantalla en la esquina inferior izquierda
        self.gui_manager.add(
            arcade.gui.UIAnchorWidget(
                anchor_x="left",
                anchor_y="bottom",
                align_x=20,
                align_y=20,
                child=v_box,
            )
        )

    def on_update(self, delta_time):
        """Lógica del juego: Se ejecuta automáticamente unas 60 veces por segundo."""

        # Obtenemos la hora actual del sistema operativo
        current_time = time.localtime()
        hours = current_time.tm_hour % 12  # Formato de 12 horas
        minutes = current_time.tm_min
        seconds = current_time.tm_sec

        # CALCULAR LOS ÁNGULOS (En Arcade, 0° es a la derecha, sumamos en sentido horario)
        # Un reloj avanza en sentido horario, por lo que restamos el ángulo para que vaya a la derecha.

        # Ángulo de los minutos: 360 grados / 60 minutos = 6 grados por minuto
        # Ajustamos el desfase: 90 grados es "las 12" en matemáticas de computación
        self.minute_angle = 90 - (minutes * 6) - (seconds * 0.1)

        # Ángulo de las horas: 360 grados / 12 horas = 30 grados por hora
        # Añadimos la fracción de la hora que representan los minutos actuales para que se mueva suavemente
        self.hour_angle = 90 - (hours * 30) - (minutes * 0.5)

        # Aplicamos los ángulos calculados a las propiedades de los Sprites
        self.minute_hand.angle = self.minute_angle
        self.hour_hand.angle = self.hour_angle

    def on_draw(self):
        """Renderizado: Dibuja todo en la pantalla en cada fotograma."""

        # Limpia la pantalla antes de volver a dibujar
        self.clear()

        # Dibujar los elementos en orden (capas de abajo hacia arriba)
        self.clock_background.draw()  # El fondo va primero
        self.hour_hand.draw()  # Las horas abajo
        self.minute_hand.draw()  # Los minutos encima

        # Dibujar la interfaz gráfica (GUI) encima de los sprites
        self.gui_manager.draw()


# --- PUNTO DE ENTRADA DEL PROGRAMA ---
def main():
    window = ClockGame()
    window.setup()
    arcade.run()


if __name__ == "__main__":
    main()
