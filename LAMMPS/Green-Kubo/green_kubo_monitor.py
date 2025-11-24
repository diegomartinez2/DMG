#!/usr/bin/env python
"""
Este script permite seguir la evolución de un cálculo Green-Kubo LAMMPS
visualizando los datos a medida que se calculan.
"""
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
from scipy.integrate import cumulative_trapezoid  # Scipy moderno, o usar cumtrapz en versiones viejas
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.animation as animation
import os
import threading
import time
import random

# =======================================================
# 1. Lógica Científica (Modelo - Green-Kubo)
# =======================================================
class GreenKuboModel:
    """
    Maneja la lectura de datos y el cálculo físico de la conductividad.
    Basado en el script 'Green-Kubo.py'.
    """
    def __init__(self):
        # Constantes Físicas y Factores (Valores por defecto)
        self.kB_eV_K = 8.617333262145e-5
        self.FACTOR_CONVERSION = 1.5975e+10 # De eV^2/(A ps K) a W/(m K)

        # Parámetros del Sistema (Modificables desde la UI)
        self.T = 80.0          # K
        self.V = 145537.0      # A^3 (Valor actualizado del input LAMMPS)
        self.dt = 0.001        # ps

        self.filename = "J0Jt_12345.dat"
        self.file_handle = None

        # Almacenamiento de datos
        self.raw_data = [] # Lista de líneas crudas para procesar incrementalmente

        # Datos procesados para gráficas
        self.time = np.array([])
        self.acf_avg = np.array([])
        self.kappa_accum = np.array([])

    def set_parameters(self, T, V, dt, filename):
        self.T = T
        self.V = V
        self.dt = dt
        self.filename = filename
        self.reset_data()

    def reset_data(self):
        self.raw_data = []
        self.time = np.array([])
        self.acf_avg = np.array([])
        self.kappa_accum = np.array([])
        if self.file_handle:
            self.file_handle.close()
            self.file_handle = None

    def read_and_process(self):
        """
        Lee nuevas líneas del archivo y actualiza los cálculos.
        Retorna True si hubo nuevos datos.
        """
        if not os.path.exists(self.filename):
            return False

        if self.file_handle is None:
            try:
                self.file_handle = open(self.filename, 'r')
                # Intentar saltar cabecera si es archivo nuevo
                # LAMMPS fix ave/correlate tiene 3 lineas de header
                for _ in range(3):
                    line = self.file_handle.readline()
                    if not line.startswith('#'): # Si no es comentario, rebobinar (caso raro)
                        self.file_handle.seek(0)
                        break
            except Exception as e:
                print(f"Error abriendo archivo: {e}")
                return False

        # Leer todas las líneas nuevas
        lines = self.file_handle.readlines()
        if not lines:
            return False

        new_steps = []
        new_acf_x = []
        new_acf_y = []
        new_acf_z = []

        for line in lines:
            # Saltar líneas de comentario o vacías
            if line.startswith('#') or not line.strip():
                continue

            try:
                parts = list(map(float, line.split()))
                # Estructura típica fix ave/correlate:
                # Col 0: TimeStep, Col 1: N_samples, Col 2: JxJx, Col 3: JyJy, Col 4: JzJz
                # Nota: El script original usaba col 1 para Jx, pero LAMMPS suele poner N_count ahí.
                # Asumiremos formato estandard LAMMPS Output:
                if len(parts) >= 5:
                    new_steps.append(parts[0])
                    new_acf_x.append(parts[2])
                    new_acf_y.append(parts[3])
                    new_acf_z.append(parts[4])
            except ValueError:
                continue

        if not new_steps:
            return False

        # Concatenar con datos existentes
        current_steps = len(self.time)
        new_time = np.array(new_steps) * self.dt

        # Promedio de las 3 componentes (Isotrópico)
        new_acf_avg = (np.array(new_acf_x) + np.array(new_acf_y) + np.array(new_acf_z)) / 3.0

        if current_steps == 0:
            self.time = new_time
            self.acf_avg = new_acf_avg
        else:
            self.time = np.concatenate((self.time, new_time))
            self.acf_avg = np.concatenate((self.acf_avg, new_acf_avg))

        # --- CÁLCULO DE GREEN-KUBO (Dinámico) ---
        # Calculamos la integral acumulativa (simulando la integral de 0 a tau)
        # Integral = cumtrapz(ACF, x=Time)

        if len(self.time) > 1:
            # Recalcular integral completa o incremental (KISS: completa es rápida para <10k puntos)
            integral_accum = cumulative_trapezoid(self.acf_avg, self.time, initial=0)

            # Prefactor GK: V / (3 * kB * T^2)
            prefactor = self.V / (3.0 * self.kB_eV_K * (self.T**2))

            # Kappa = Prefactor * Integral * Conversión
            self.kappa_accum = prefactor * integral_accum * self.FACTOR_CONVERSION

        return True

# =======================================================
# 2. Simulador de Datos LAMMPS (Para pruebas)
# =======================================================
class FakeLAMMPS(threading.Thread):
    def __init__(self, filename="J0Jt_12345.dat"):
        super().__init__()
        self.filename = filename
        self.running = False
        self.daemon = True # Morir si el programa principal muere

    def run(self):
        self.running = True
        with open(self.filename, 'w') as f:
            f.write("# TimeStep N-count c_myFlux[1] c_myFlux[2] c_myFlux[3]\n")
            f.write("# Simulated Data for Testing\n")
            f.write("# \n")

        step = 0
        dt_step = 10
        # Generar una curva de decaimiento exponencial con ruido (ACF típica)
        tau_decay = 200 # pasos

        while self.running:
            time.sleep(0.1) # Velocidad de escritura simulada

            # Modelo simple de ACF: A * exp(-t/tau) + ruido
            decay = np.exp(-step / tau_decay)

            # Generar 3 componentes con ruido correlacionado
            jx = 5.0 * decay + random.uniform(-0.5, 0.5) * decay
            jy = 5.0 * decay + random.uniform(-0.5, 0.5) * decay
            jz = 5.0 * decay + random.uniform(-0.5, 0.5) * decay

            n_count = 1000

            with open(self.filename, 'a') as f:
                f.write(f"{step} {n_count} {jx:.6f} {jy:.6f} {jz:.6f}\n")

            step += dt_step

            if step > 5000: # Reiniciar o detener después de un tiempo
                 # self.running = False
                 pass

    def stop(self):
        self.running = False

# =======================================================
# 3. Interfaz Gráfica (Vista)
# =======================================================
class GKMonitorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Monitor de Conductividad Térmica (Green-Kubo)")
        self.root.geometry("1100x800")

        self.model = GreenKuboModel()
        self.simulator = None
        self.is_monitoring = False

        self._setup_ui()
        self._init_plots()

        # Animación para refresco
        self.ani = animation.FuncAnimation(self.fig, self.update_plots, interval=1000, cache_frame_data=False)

    def _setup_ui(self):
        # --- Panel de Control (Izquierda) ---
        control_frame = tk.Frame(self.root, bg="#e1e1e1", padx=10, pady=10, width=250)
        control_frame.pack(side=tk.LEFT, fill=tk.Y)
        control_frame.pack_propagate(False) # Mantener ancho fijo

        tk.Label(control_frame, text="Configuración", font=("Arial", 14, "bold"), bg="#e1e1e1").pack(pady=10)

        # Entradas de Parámetros
        self.ent_T = self._add_input(control_frame, "Temperatura (K):", "80.0")
        self.ent_V = self._add_input(control_frame, "Volumen (A^3):", "145537.0")
        self.ent_dt = self._add_input(control_frame, "Timestep (ps):", "0.001")
        self.ent_file = self._add_input(control_frame, "Archivo ACF:", "J0Jt_12345.dat")

        # Botón Selector de Archivo
        tk.Button(control_frame, text="Buscar Archivo", command=self.browse_file).pack(fill=tk.X, pady=5)

        ttk.Separator(control_frame, orient='horizontal').pack(fill='x', pady=20)

        # Botones de Acción
        self.btn_monitor = tk.Button(control_frame, text="Iniciar Monitorización",
                                     bg="green", fg="white", font=("Arial", 11, "bold"),
                                     command=self.toggle_monitoring)
        self.btn_monitor.pack(fill=tk.X, pady=10)

        ttk.Separator(control_frame, orient='horizontal').pack(fill='x', pady=20)

        tk.Label(control_frame, text="Pruebas / Demo", font=("Arial", 12, "bold"), bg="#e1e1e1").pack(pady=5)
        self.btn_sim = tk.Button(control_frame, text="Generar Datos Falsos",
                                 bg="#555", fg="white", command=self.toggle_simulation)
        self.btn_sim.pack(fill=tk.X, pady=5)

        # Panel de Resultados Numéricos
        self.lbl_kappa_curr = tk.Label(control_frame, text="Kappa Actual:\n--- W/mK",
                                       font=("Arial", 16, "bold"), bg="#e1e1e1", fg="#0000aa")
        self.lbl_kappa_curr.pack(side=tk.BOTTOM, pady=20)

    def _add_input(self, parent, label, default):
        tk.Label(parent, text=label, bg="#e1e1e1", anchor="w").pack(fill=tk.X)
        entry = tk.Entry(parent)
        entry.insert(0, default)
        entry.pack(fill=tk.X, pady=(0, 10))
        return entry

    def _init_plots(self):
        # --- Área de Gráficos (Derecha) ---
        plot_frame = tk.Frame(self.root, bg="white")
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        self.fig = Figure(figsize=(8, 8), dpi=100)

        # Subplot 1: Autocorrelación (ACF)
        self.ax1 = self.fig.add_subplot(211)
        self.ax1.set_title("Autocorrelación del Flujo de Calor (ACF)")
        self.ax1.set_ylabel("ACF (Promedio)")
        self.ax1.axhline(0, color='gray', linestyle='--', linewidth=0.8)
        self.ax1.grid(True, alpha=0.5)

        # Subplot 2: Integral (Kappa)
        self.ax2 = self.fig.add_subplot(212)
        self.ax2.set_title("Conductividad Térmica Convergente (Integral)")
        self.ax2.set_xlabel("Tiempo de Correlación (ps)")
        self.ax2.set_ylabel("Kappa (W/mK)")
        self.ax2.grid(True, alpha=0.5)

        self.fig.tight_layout(pad=3.0)

        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def browse_file(self):
        f = filedialog.askopenfilename(filetypes=[("Data files", "*.dat"), ("All files", "*.*")])
        if f:
            self.ent_file.delete(0, tk.END)
            self.ent_file.insert(0, f)

    def toggle_simulation(self):
        if self.simulator and self.simulator.is_alive():
            self.simulator.stop()
            self.btn_sim.config(text="Generar Datos Falsos", bg="#555")
        else:
            # Usar el nombre de archivo actual
            fname = self.ent_file.get()
            self.simulator = FakeLAMMPS(fname)
            self.simulator.start()
            self.btn_sim.config(text="Detener Generación", bg="orange")
            # Auto-iniciar monitoreo si no está activo
            if not self.is_monitoring:
                self.toggle_monitoring()

    def toggle_monitoring(self):
        if self.is_monitoring:
            self.is_monitoring = False
            self.btn_monitor.config(text="Iniciar Monitorización", bg="green")
        else:
            # Aplicar parámetros
            try:
                T = float(self.ent_T.get())
                V = float(self.ent_V.get())
                dt = float(self.ent_dt.get())
                fname = self.ent_file.get()

                self.model.set_parameters(T, V, dt, fname)
                self.is_monitoring = True
                self.btn_monitor.config(text="Detener Monitorización", bg="red")
            except ValueError:
                messagebox.showerror("Error", "Por favor revisa que los parámetros numéricos sean válidos.")

    def update_plots(self, frame):
        if not self.is_monitoring:
            return

        # 1. Leer datos
        has_new_data = self.model.read_and_process()

        if has_new_data and len(self.model.time) > 0:
            # 2. Actualizar Gráficos
            self.ax1.clear()
            self.ax2.clear()

            # ACF Plot
            self.ax1.plot(self.model.time, self.model.acf_avg, color='#444444', linewidth=1, label='ACF Promedio')
            self.ax1.set_title("Autocorrelación del Flujo de Calor")
            self.ax1.set_ylabel("Unidades LAMMPS")
            self.ax1.axhline(0, color='red', linestyle='--', alpha=0.5)
            self.ax1.legend(loc='upper right')
            self.ax1.grid(True, alpha=0.3)

            # Kappa Plot
            # Graficamos la integral acumulativa. Se busca el "plateau" (meseta).
            if len(self.model.kappa_accum) > 0:
                # Ajustamos x para que coincida con kappa (cumtrapz devuelve N-1 elementos)
                t_plot = self.model.time[1:]
                k_plot = self.model.kappa_accum

                self.ax2.plot(t_plot, k_plot, color='blue', linewidth=2, label='Kappa (Integral)')
                self.ax2.set_title("Convergencia de Conductividad Térmica")
                self.ax2.set_xlabel("Tiempo de Correlación (ps)")
                self.ax2.set_ylabel("Kappa (W/m K)")
                self.ax2.grid(True, alpha=0.3)
                self.ax2.legend(loc='lower right')

                # Actualizar etiqueta numérica con el último valor
                last_val = k_plot[-1]
                self.lbl_kappa_curr.config(text=f"Kappa Actual:\n{last_val:.3f} W/mK")

            self.canvas.draw()

    def on_close(self):
        self.is_monitoring = False
        if self.simulator:
            self.simulator.stop()
        self.root.destroy()

if __name__ == "__main__":
    root = tk.Tk()
    app = GKMonitorApp(root)
    root.protocol("WM_DELETE_WINDOW", app.on_close)
    root.mainloop()
