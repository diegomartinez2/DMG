import tkinter as tk
from tkinter import ttk
import random
import time

class HabiJumpInterface:
    def __init__(self, root):
        self.root = root
        self.root.title("HABI OS - Módulo de Salto de Fase v3.0")
        self.root.geometry("600x700")
        self.root.configure(bg="#0a0a0c")

        # Configuración de estilos
        self.style = ttk.Style()
        self.style.theme_use('clam')

        # Colores
        self.accent_color = "#00ffcc"  # Cian Habi
        self.bg_color = "#0a0a0c"
        self.warn_color = "#ff3366"

        self.create_widgets()
        self.update_simulation()

    def create_widgets(self):
        # Header
        header_frame = tk.Frame(self.root, bg=self.bg_color, pady=20)
        header_frame.pack(fill="x")

        tk.Label(header_frame, text="HABI INTERFACE SYSTEM", font=("Courier", 10), fg=self.accent_color, bg=self.bg_color).pack()
        tk.Label(header_frame, text="TRANSCODER DE FASE PSICOMÉTRICO", font=("Orbitron", 16, "bold"), fg="white", bg=self.bg_color).pack()

        # Monitor de Estado de Firma
        self.signature_frame = tk.LabelFrame(self.root, text=" FIRMA PSICOMÉTRICA ", font=("Courier", 9), fg=self.accent_color, bg=self.bg_color, padx=15, pady=15)
        self.signature_frame.pack(fill="x", padx=20, pady=10)

        self.sig_status = tk.Label(self.signature_frame, text="SINCRONIZANDO...", font=("Courier", 12, "bold"), fg="#ffff00", bg=self.bg_color)
        self.sig_status.pack(side="left")

        self.sig_bar = ttk.Progressbar(self.signature_frame, length=200, mode='determinate')
        self.sig_bar.pack(side="right")

        # Radar de Solapamiento (Canvas)
        radar_frame = tk.Frame(self.root, bg=self.bg_color)
        radar_frame.pack(pady=10)

        tk.Label(radar_frame, text="ESCÁNER DE SOLAPAMIENTO MOLECULAR", font=("Courier", 9), fg=self.accent_color, bg=self.bg_color).pack()
        self.canvas = tk.Canvas(radar_frame, width=300, height=300, bg="#001111", highlightthickness=1, highlightbackground=self.accent_color)
        self.canvas.pack()
        self.draw_radar_base()

        # Condensadores Inerciales
        condenser_frame = tk.Frame(self.root, bg=self.bg_color, padx=20)
        condenser_frame.pack(fill="x", pady=20)

        tk.Label(condenser_frame, text="BANCO DE CONDENSADORES (ENERGÍA CINÉTICA)", font=("Courier", 9), fg=self.accent_color, bg=self.bg_color).pack(anchor="w")
        self.energy_val = tk.IntVar(value=65)
        self.energy_bar = ttk.Progressbar(condenser_frame, variable=self.energy_val, maximum=100)
        self.energy_bar.pack(fill="x", pady=5)
        self.energy_label = tk.Label(condenser_frame, text="65% CAPACIDAD", font=("Courier", 10), fg="white", bg=self.bg_color)
        self.energy_label.pack(anchor="e")

        # Botón de Salto
        self.jump_btn = tk.Button(self.root, text="INICIAR SALTO DE FASE", font=("Orbitron", 14, "bold"),
                                  bg=self.accent_color, fg=self.bg_color, activebackground="white",
                                  command=self.execute_jump, cursor="hand2", bd=0, padx=20, pady=10)
        self.jump_btn.pack(pady=20)

        # Consola de Mensajes
        self.console = tk.Text(self.root, height=4, bg="black", fg="#00ff00", font=("Courier", 8), bd=0, padx=10, pady=10)
        self.console.pack(fill="x", padx=20, pady=10)
        self.log("Sistemas en línea. Esperando comando del usuario.")

    def draw_radar_base(self):
        # Dibuja círculos de radar
        for i in range(1, 4):
            r = i * 50
            self.canvas.create_oval(150-r, 150-r, 150+r, 150+r, outline="#004444")
        self.canvas.create_line(150, 0, 150, 300, fill="#004444")
        self.canvas.create_line(0, 150, 300, 150, fill="#004444")
        self.user_dot = self.canvas.create_oval(145, 145, 155, 155, fill=self.accent_color, outline="white")

    def log(self, message):
        self.console.insert("1.0", f"> {message}\n")

    def update_simulation(self):
        # Simular fluctuación de firma
        val = random.randint(95, 100)
        self.sig_bar['value'] = val
        if val > 98:
            self.sig_status.config(text="FIRMA ESTABLE", fg=self.accent_color)
        else:
            self.sig_status.config(text="RE-SINCRONIZANDO", fg="#ffff00")

        # Simular interferencias en radar
        self.canvas.delete("obs")
        for _ in range(3):
            x = random.randint(50, 250)
            y = random.randint(50, 250)
            self.canvas.create_oval(x, y, x+5, y+5, fill=self.warn_color, tags="obs")

        self.root.after(1000, self.update_simulation)

    def execute_jump(self):
        self.jump_btn.config(state="disabled", text="DESFASANDO...")
        self.log("Iniciando secuencia de colapso de onda...")
        self.log("Validando coordenadas con Habi...")

        # Efecto de carga
        def complete():
            self.log("¡SALTO EXITOSO! Relocalización completada.")
            self.energy_val.set(self.energy_val.get() + random.randint(-10, 10))
            self.energy_label.config(text=f"{self.energy_val.get()}% CAPACIDAD")
            self.jump_btn.config(state="normal", text="INICIAR SALTO DE FASE")

        self.root.after(2000, complete)

if __name__ == "__main__":
    root = tk.Tk()
    app = HabiJumpInterface(root)
    root.mainloop()
