import pygame
import sys

# Configuración inicial
WIDTH, HEIGHT = 1024, 768
FPS = 60

# Colores (Paleta industrial/retro-futurista)
COLOR_BG = (20, 25, 30)
COLOR_PANEL = (45, 55, 65)
COLOR_TEXT = (0, 255, 65)  # Verde fósforo
COLOR_BTN = (70, 80, 90)
COLOR_BTN_HOVER = (100, 110, 120)
COLOR_DANGER = (255, 50, 50)

class Button:
    def __init__(self, x, y, w, h, text, callback):
        self.rect = pygame.Rect(x, y, w, h)
        self.text = text
        self.callback = callback
        self.hovered = False

    def draw(self, screen, font):
        color = COLOR_BTN_HOVER if self.hovered else COLOR_BTN
        pygame.draw.rect(screen, color, self.rect, border_radius=5)
        pygame.draw.rect(screen, COLOR_TEXT, self.rect, 2, border_radius=5)

        text_surf = font.render(self.text, True, COLOR_TEXT)
        text_rect = text_surf.get_rect(center=self.rect.center)
        screen.blit(text_surf, text_rect)

    def check_hover(self, pos):
        self.hovered = self.rect.collidepoint(pos)

    def handle_event(self, event):
        if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
            if self.hovered:
                self.callback()

class App:
    def __init__(self):
        pygame.init()
        self.screen = pygame.display.set_caption("NSS-Console: Sistema de Gestión de Nave")
        self.screen = pygame.display.set_mode((WIDTH, HEIGHT))
        self.clock = pygame.time.Clock()
        self.font_main = pygame.font.SysFont("Courier New", 18, bold=True)
        self.font_ui = pygame.font.SysFont("Courier New", 14)

        # Estado del juego
        self.current_station = "REACTOR" # REACTOR, SUPPORT, DISTRO
        self.reactor_temp = 45.0
        self.energy_output = 100
        self.coolant_active = False

        # Crear botones de navegación
        self.nav_buttons = [
            Button(20, 20, 150, 40, "1. REACTOR", lambda: self.set_station("REACTOR")),
            Button(180, 20, 150, 40, "2. SOPORTE VITAL", lambda: self.set_station("SUPPORT")),
            Button(340, 20, 150, 40, "3. DISTRIBUCION", lambda: self.set_station("DISTRO")),
        ]

        # Botones específicos
        self.control_buttons = {
            "REACTOR": [
                Button(400, 500, 200, 50, "BOMBA REFRIG.", self.toggle_coolant)
            ]
        }

    def set_station(self, name):
        self.current_station = name

    def toggle_coolant(self):
        self.coolant_active = not self.coolant_active

    def update_logic(self):
        # Lógica del simulador
        if self.coolant_active:
            self.reactor_temp -= 0.1
        else:
            self.reactor_temp += 0.05

        # Limites
        if self.reactor_temp < 20: self.reactor_temp = 20
        if self.reactor_temp > 200: self.reactor_temp = 200

    def draw_station_reactor(self):
        # Dibujar Panel Central
        pygame.draw.rect(self.screen, COLOR_PANEL, (50, 100, 924, 600))

        # Indicador de Temperatura
        label = self.font_main.render(f"TEMPERATURA DEL NUCLEO: {self.reactor_temp:.2f} C", True, COLOR_TEXT)
        self.screen.blit(label, (100, 150))

        # Barra de peligro
        bar_width = 400
        fill_width = (self.reactor_temp / 200) * bar_width
        pygame.draw.rect(self.screen, (30, 30, 30), (100, 180, bar_width, 30))
        color = COLOR_DANGER if self.reactor_temp > 150 else COLOR_TEXT
        pygame.draw.rect(self.screen, color, (100, 180, fill_width, 30))

        # Estado bomba
        status_color = (0, 255, 0) if self.coolant_active else (100, 0, 0)
        status_text = "ACTIVA" if self.coolant_active else "APAGADA"
        st_label = self.font_ui.render(f"BOMBA REFRIGERANTE: {status_text}", True, status_color)
        self.screen.blit(st_label, (100, 220))

    def run(self):
        while True:
            mouse_pos = pygame.mouse.get_pos()

            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    pygame.quit()
                    sys.exit()

                for btn in self.nav_buttons:
                    btn.handle_event(event)

                if self.current_station in self.control_buttons:
                    for btn in self.control_buttons[self.current_station]:
                        btn.handle_event(event)

            # Actualizar hovers
            for btn in self.nav_buttons: btn.check_hover(mouse_pos)
            if self.current_station in self.control_buttons:
                for btn in self.control_buttons[self.current_station]: btn.check_hover(mouse_pos)

            self.update_logic()

            # Renderizado
            self.screen.fill(COLOR_BG)

            # Cabecera
            pygame.draw.rect(self.screen, (10, 10, 15), (0, 0, WIDTH, 80))
            for btn in self.nav_buttons:
                btn.draw(self.screen, self.font_ui)

            # Dibujar estación actual
            if self.current_station == "REACTOR":
                self.draw_station_reactor()
                for btn in self.control_buttons["REACTOR"]:
                    btn.draw(self.screen, self.font_main)
            else:
                msg = self.font_main.render(f"ESTACION: {self.current_station} - EN DESARROLLO", True, COLOR_TEXT)
                self.screen.blit(msg, (WIDTH//2 - 200, HEIGHT//2))

            pygame.display.flip()
            self.clock.tick(FPS)

if __name__ == "__main__":
    App().run()
