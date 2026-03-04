import arcade
import arcade.gui

# Constantes de pantalla
SCREEN_WIDTH = 1024
SCREEN_HEIGHT = 768
SCREEN_TITLE = "Sistema Operativo de Nave Espacial (Arcade Edition)"

# Colores Estilo Terminal
COLOR_CRTR_BG = (10, 15, 10)
COLOR_UI_BORDER = (0, 255, 65)

class StationView(arcade.View):
    """Clase base para todas las estaciones de control"""
    def __init__(self):
        super().__init__()
        self.manager = arcade.gui.UIManager()

    def on_show_view(self):
        self.manager.enable()
        arcade.set_background_color(COLOR_CRTR_BG)

    def on_hide_view(self):
        self.manager.disable()

    def draw_header(self, title):
        # Dibujar marco de la pantalla
        arcade.draw_rectangle_outline(
            SCREEN_WIDTH // 2, SCREEN_HEIGHT // 2,
            SCREEN_WIDTH - 20, SCREEN_HEIGHT - 20,
            COLOR_UI_BORDER, 2
        )
        arcade.draw_text(title, SCREEN_WIDTH // 2, SCREEN_HEIGHT - 50,
                         arcade.color.GREEN, font_size=24, anchor_x="center")

class ReactorView(StationView):
    """Panel de Control del Reactor"""
    def __init__(self):
        super().__init__()
        self.temp = 50.0

        # UI Manager para botones - space_between gestiona el espacio entre elementos
        self.v_box = arcade.gui.UIBoxLayout(space_between=20)

        # Botón para cambiar a Soporte Vital
        btn_support = arcade.gui.UIFlatButton(text="Ir a Soporte Vital", width=200)
        self.v_box.add(btn_support)

        # Evento de click (usando decorador moderno)
        @btn_support.event("on_click")
        def on_click_flatbutton(event):
            self.window.show_view(LifeSupportView())

        # Centrar el layout en la parte inferior
        self.manager.add(
            arcade.gui.UIAnchorWidget(
                anchor_x="center_x",
                anchor_y="bottom",
                align_y=50,
                child=self.v_box)
        )

    def on_update(self, delta_time):
        self.temp += delta_time * 2

    def on_draw(self):
        self.clear()
        self.draw_header("ESTACION 01: NUCLEO DEL REACTOR")

        status_color = arcade.color.GREEN if self.temp < 100 else arcade.color.RED
        arcade.draw_text(f"TEMPERATURA: {self.temp:.1f}°K", 100, 500,
                         status_color, font_size=20)

        # Gráfica básica
        arcade.draw_line(100, 200, 900, 200, COLOR_UI_BORDER, 2)
        arcade.draw_line(100, 200, 100, 450, COLOR_UI_BORDER, 2)

        self.manager.draw()

class LifeSupportView(StationView):
    """Panel de Soporte Vital"""
    def __init__(self):
        super().__init__()

        self.v_box = arcade.gui.UIBoxLayout()
        btn_reactor = arcade.gui.UIFlatButton(text="Volver a Reactor", width=200)
        self.v_box.add(btn_reactor)

        @btn_reactor.event("on_click")
        def on_click_flatbutton(event):
            self.window.show_view(ReactorView())

        self.manager.add(
            arcade.gui.UIAnchorWidget(anchor_x="center_x", anchor_y="bottom", align_y=50, child=self.v_box)
        )

    def on_draw(self):
        self.clear()
        self.draw_header("ESTACION 02: SOPORTE VITAL")

        arcade.draw_text("OXIGENO: 98%", 100, 500, arcade.color.AERO_BLUE, font_size=20)
        arcade.draw_text("PRESION: 1.0 ATM", 100, 460, arcade.color.AERO_BLUE, font_size=20)

        self.manager.draw()

def main():
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE, resizable=True)
    start_view = ReactorView()
    window.show_view(start_view)
    arcade.run()

if __name__ == "__main__":
    main()
