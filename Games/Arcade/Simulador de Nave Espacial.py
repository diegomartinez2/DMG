import arcade
import arcade.gui

# Constantes de pantalla
SCREEN_WIDTH = 1024
SCREEN_HEIGHT = 768
SCREEN_TITLE = "Sistema Operativo de Nave Espacial"

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
        # Usamos arcade.rect.XYWH o arcade.rect.LRBT para crear el rectángulo correctamente
        # LRBT: Left, Right, Bottom, Top
        border_rect = arcade.rect.LRBT(
            left=10,
            right=self.window.width - 10,
            bottom=10,
            top=self.window.height - 10
        )

        # Dibujamos el contorno
        arcade.draw_rect_outline(
            rect=border_rect,
            color=COLOR_UI_BORDER,
            border_width=2
        )

        arcade.draw_text(
            title,
            self.window.width // 2,
            self.window.height - 50,
            arcade.color.GREEN,
            font_size=24,
            anchor_x="center"
        )

class ReactorView(StationView):
    def __init__(self):
        super().__init__()
        self.temp = 50.0

        # Layout principal
        self.anchor = arcade.gui.UIAnchorLayout()
        self.v_box = arcade.gui.UIBoxLayout(space_between=20)

        # Botón para cambiar de vista
        btn_support = arcade.gui.UIFlatButton(text="IR A SOPORTE VITAL", width=250)
        self.v_box.add(btn_support)

        @btn_support.event("on_click")
        def on_click_btn(event):
            self.window.show_view(LifeSupportView())

        self.anchor.add(
            child=self.v_box,
            anchor_x="center_x",
            anchor_y="bottom",
            align_y=50
        )
        self.manager.add(self.anchor)

    def on_update(self, delta_time):
        # Simulación simple de temperatura
        self.temp += delta_time * 1.5

    def on_draw(self):
        self.clear()
        self.draw_header("ESTACION 01: NUCLEO DEL REACTOR")

        status_color = arcade.color.GREEN if self.temp < 80 else arcade.color.ORANGE_PEEL
        if self.temp > 110: status_color = arcade.color.RED

        arcade.draw_text(f"ESTADO DEL NUCLEO: OPERATIVO", 100, 550, arcade.color.GREEN, font_size=18)
        arcade.draw_text(f"TEMPERATURA: {self.temp:.1f}°K", 100, 510, status_color, font_size=22, bold=True)

        # Líneas decorativas
        arcade.draw_line(100, 480, 400, 480, COLOR_UI_BORDER, 2)

        self.manager.draw()

class LifeSupportView(StationView):
    def __init__(self):
        super().__init__()
        self.anchor = arcade.gui.UIAnchorLayout()
        self.v_box = arcade.gui.UIBoxLayout()

        btn_reactor = arcade.gui.UIFlatButton(text="VOLVER AL REACTOR", width=250)
        self.v_box.add(btn_reactor)

        @btn_reactor.event("on_click")
        def on_click_btn(event):
            self.window.show_view(ReactorView())

        self.anchor.add(
            child=self.v_box,
            anchor_x="center_x",
            anchor_y="bottom",
            align_y=50
        )
        self.manager.add(self.anchor)

    def on_draw(self):
        self.clear()
        self.draw_header("ESTACION 02: SOPORTE VITAL")

        arcade.draw_text("NIVELES DE OXIGENO: 98.4%", 100, 550, arcade.color.AERO_BLUE, font_size=20)
        arcade.draw_text("PRESION CABINA: 1.02 ATM", 100, 510, arcade.color.AERO_BLUE, font_size=20)
        arcade.draw_text("FILTROS CO2: NOMINAL", 100, 470, arcade.color.GREEN, font_size=18)

        self.manager.draw()

def main():
    window = arcade.Window(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE, resizable=True)
    # Iniciamos con la vista del Reactor
    window.show_view(ReactorView())
    arcade.run()

if __name__ == "__main__":
    main()
