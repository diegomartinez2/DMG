import arcade
import arcade.gui

# --- CONFIGURACIÓN ---
SCREEN_WIDTH = 1024
SCREEN_HEIGHT = 768
SCREEN_TITLE = "Sistema Operativo de Nave Espacial - v0.2"

# Estética Terminal
COLOR_CRTR_BG = (10, 15, 10)
COLOR_UI_BORDER = (0, 255, 65)

class GameWindow(arcade.Window):
    """
    Clase principal que mantiene el estado global del juego.
    Los datos aquí definidos persisten durante toda la sesión.
    """
    def __init__(self):
        super().__init__(SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_TITLE, resizable=True)

        # --- ESTADO GLOBAL DE LA NAVE ---
        self.ship_stats = {
            "reactor_temp": 50.0,
            "oxygen_level": 98.4,
            "hull_integrity": 100,
            "power_output": 75.0
        }

    def on_update(self, delta_time):
        # El reactor siempre se calienta, sin importar en qué vista estemos
        self.ship_stats["reactor_temp"] += delta_time * 0.8

        # Simulación de consumo de oxígeno
        if self.ship_stats["oxygen_level"] > 0:
            self.ship_stats["oxygen_level"] -= delta_time * 0.05

class StationView(arcade.View):
    """Clase base con utilidades de dibujo comunes"""
    def __init__(self):
        super().__init__()
        self.manager = arcade.gui.UIManager()

    def on_show_view(self):
        self.manager.enable()
        arcade.set_background_color(COLOR_CRTR_BG)

    def on_hide_view(self):
        self.manager.disable()

    def draw_ui_frame(self, title):
        # Marco principal persistente
        border_rect = arcade.rect.LRBT(
            left=15,
            right=self.window.width - 15,
            bottom=15,
            top=self.window.height - 15
        )
        arcade.draw_rect_outline(rect=border_rect, color=COLOR_UI_BORDER, border_width=2)

        # Título de la Estación
        arcade.draw_text(
            title,
            self.window.width // 2,
            self.window.height - 60,
            arcade.color.GREEN,
            font_size=28,
            anchor_x="center",
            font_name="Courier New"
        )

        # Barra de estado inferior rápida (Mini-HUD)
        stats = self.window.ship_stats
        hud_text = f"TEMP: {stats['reactor_temp']:.1f}°K | O2: {stats['oxygen_level']:.1f}% | CASCO: {stats['hull_integrity']}%"
        arcade.draw_text(hud_text, 30, 30, COLOR_UI_BORDER, font_size=12)

class ReactorView(StationView):
    def __init__(self):
        super().__init__()

        self.anchor = arcade.gui.UIAnchorLayout()
        self.v_box = arcade.gui.UIBoxLayout(space_between=20)

        # Navegación
        btn_support = arcade.gui.UIFlatButton(text="SOPORTE VITAL", width=200)
        @btn_support.event("on_click")
        def on_click_btn(event):
            self.window.show_view(LifeSupportView())

        self.v_box.add(btn_support)
        self.anchor.add(child=self.v_box, anchor_x="right", anchor_y="bottom", align_x=-50, align_y=50)
        self.manager.add(self.anchor)

    def on_draw(self):
        self.clear()
        self.draw_ui_frame("NUCLEO DEL REACTOR")

        temp = self.window.ship_stats["reactor_temp"]

        # Color dinámico según peligro
        color = arcade.color.GREEN
        if temp > 100: color = arcade.color.RED
        elif temp > 80: color = arcade.color.ORANGE

        # Dibujo de la "barra de calor"
        arcade.draw_text("NIVEL TÉRMICO:", 100, 500, arcade.color.WHITE, font_size=16)
        arcade.draw_lrbt_rectangle_filled(100, 100 + (temp * 2), 450, 480, color)
        arcade.draw_lrbt_rectangle_outline(100, 400, 450, 480, arcade.color.GRAY, 2)

        arcade.draw_text(f"{temp:.2f} °K", 420, 455, color, font_size=20, bold=True)

        if temp > 100:
            arcade.draw_text("!!! ALERTA DE SOBRECALENTAMIENTO !!!", self.window.width//2, 300,
                             arcade.color.RED, font_size=20, anchor_x="center")

        self.manager.draw()

class LifeSupportView(StationView):
    def __init__(self):
        super().__init__()
        self.anchor = arcade.gui.UIAnchorLayout()

        btn_reactor = arcade.gui.UIFlatButton(text="REGRESAR AL REACTOR", width=200)
        @btn_reactor.event("on_click")
        def on_click_btn(event):
            self.window.show_view(ReactorView())

        self.anchor.add(child=btn_reactor, anchor_x="right", anchor_y="bottom", align_x=-50, align_y=50)
        self.manager.add(self.anchor)

    def on_draw(self):
        self.clear()
        self.draw_ui_frame("SISTEMAS DE SOPORTE VITAL")

        o2 = self.window.ship_stats["oxygen_level"]
        arcade.draw_text(f"SUMINISTRO DE OXÍGENO: {o2:.2f}%", 100, 500, arcade.color.AERO_BLUE, font_size=22)

        self.manager.draw()

def main():
    # Creamos la ventana personalizada que guarda el estado
    window = GameWindow()
    window.show_view(ReactorView())
    arcade.run()

if __name__ == "__main__":
    main()
