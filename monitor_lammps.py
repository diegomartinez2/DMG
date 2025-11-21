import tkinter as tk
from tkinter import ttk
import time
import random
import threading
import os

# --- Importaciones de Matplotlib ---
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.animation as animation

# --- SIMULADOR DE LAMMPS (Para pruebas) ---
class LammpsSimulator(threading.Thread):
    """
    Simula un proceso de LAMMPS escribiendo datos en un archivo.
    Se ejecuta en un hilo separado para no bloquear la UI.
    """
    def __init__(self, filename="lammps_log.txt"):
        super().__init__()
        self.filename = filename
        self.running = True

    def run(self):
        # Limpiar/Crear archivo
        with open(self.filename, "w") as f:
            f.write("Step Temp PotEng KinEng TotEng Press Volume\n")

        step = 0
        temp = 300.0

        while self.running:
            time.sleep(0.5) # Simula el tiempo de cálculo de LAMMPS

            # Generar datos con un poco de ruido (random walk)
            step += 100
            temp += random.uniform(-5, 5)
            pot_eng = -25000 + random.uniform(-10, 10)
            kin_eng = temp * 1.5 * 8.314  # Aproximación simple (3/2 N kT)
            tot_eng = pot_eng + kin_eng
            press = 1.0 + random.uniform(-0.1, 0.1)
            vol = 1000.0

            line = f"{step} {temp:.2f} {pot_eng:.2f} {kin_eng:.2f} {tot_eng:.2f} {press:.2f} {vol:.2f}\n"

            with open(self.filename, "a") as f:
                f.write(line)
                f.flush() # Forzar escritura a disco

    def stop(self):
        self.running = False


# --- LECTOR DE LOGS (Modelo - SRP) ---
class LogReader:
    """
    Se encarga únicamente de leer las nuevas líneas del archivo.
    """
    def __init__(self, filename):
        self.filename = filename
        self.file_handle = None
        self.steps = []
        self.temps = []
        self.tot_engs = []

    def read_new_data(self):
        """Lee solo las líneas nuevas añadidas al archivo."""
        if not os.path.exists(self.filename):
            return

        if self.file_handle is None:
            try:
                self.file_handle = open(self.filename, "r")
                # Saltamos la cabecera si es la primera vez
                self.file_handle.readline()
            except IOError:
                return

        # Leer todas las líneas nuevas disponibles
        lines = self.file_handle.readlines()

        for line in lines:
            parts = line.split()
            if len(parts) < 5: continue # Línea incompleta o vacía

            try:
                # Asumimos formato: Step Temp PotEng KinEng TotEng ...
                s = int(parts[0])
                t = float(parts[1])
                e = float(parts[4])

                self.steps.append(s)
                self.temps.append(t)
                self.tot_engs.append(e)
            except ValueError:
                continue # Ignorar líneas de texto no numérico

        return self.steps, self.temps, self.tot_engs

    def close(self):
        if self.file_handle:
            self.file_handle.close()


# --- VISUALIZADOR (Vista/Controlador) ---
class LiveGraphApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Monitor LAMMPS en Vivo")
        self.root.geometry("1000x600")

        # 1. Iniciar simulador y lector
        self.log_file = "lammps_log.txt"
        self.simulator = LammpsSimulator(self.log_file)
        self.reader = LogReader(self.log_file)

        # 2. UI Layout
        self._init_ui()

        # 3. Configuración de la Animación de Matplotlib
        # FuncAnimation llama a self.update_graph cada 1000ms (1s)
        self.ani = animation.FuncAnimation(self.fig, self.update_graph, interval=1000, cache_frame_data=False)

    def _init_ui(self):
        # Botones de Control
        frame_controls = tk.Frame(self.root, bg="#ddd")
        frame_controls.pack(side=tk.TOP, fill=tk.X)

        tk.Button(frame_controls, text="Iniciar Simulación (Demo)",
                  command=self.start_simulation, bg="green", fg="white").pack(side=tk.LEFT, padx=10, pady=5)

        tk.Button(frame_controls, text="Detener Simulación",
                  command=self.stop_simulation, bg="red", fg="white").pack(side=tk.LEFT, padx=10, pady=5)

        # Área de Gráficas
        self.fig = Figure(figsize=(10, 6), dpi=100)

        # Subplot 1: Temperatura
        self.ax1 = self.fig.add_subplot(211) # 2 filas, 1 col, pos 1
        self.ax1.set_title("Temperatura del Sistema")
        self.ax1.set_ylabel("Temp (K)")
        self.ax1.grid(True)

        # Subplot 2: Energía Total
        self.ax2 = self.fig.add_subplot(212) # 2 filas, 1 col, pos 2
        self.ax2.set_title("Energía Total")
        self.ax2.set_xlabel("Paso de Tiempo (Step)")
        self.ax2.set_ylabel("Energía (eV)")
        self.ax2.grid(True)

        self.fig.tight_layout() # Ajustar márgenes

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def start_simulation(self):
        if not self.simulator.is_alive():
            self.simulator = LammpsSimulator(self.log_file)
            self.simulator.start()

    def stop_simulation(self):
        if self.simulator.is_alive():
            self.simulator.stop()

    def update_graph(self, frame):
        """Esta función es llamada automáticamente por FuncAnimation."""

        # 1. Leer nuevos datos del archivo
        data = self.reader.read_new_data()
        if not data: return # No hay archivo o datos nuevos

        steps, temps, engs = data

        # 2. Limpiar y Graficar (KISS: redibujar es más fácil que actualizar líneas)
        self.ax1.clear()
        self.ax2.clear()

        # Re-aplicar estilos (clear borra todo)
        self.ax1.set_title("Temperatura")
        self.ax1.set_ylabel("K")
        self.ax1.grid(True)

        self.ax2.set_title("Energía Total")
        self.ax2.set_ylabel("eV")
        self.ax2.grid(True)

        # Graficar los últimos 500 puntos para que no se sature
        limit = -500
        self.ax1.plot(steps[limit:], temps[limit:], 'b-')
        self.ax2.plot(steps[limit:], engs[limit:], 'r-')

    def on_closing(self):
        self.stop_simulation()
        self.reader.close()
        self.root.destroy()

if __name__ == "__main__":
    root = tk.Tk()
    app = LiveGraphApp(root)
    root.protocol("WM_DELETE_WINDOW", app.on_closing) # Manejar cierre limpio
    root.mainloop()
