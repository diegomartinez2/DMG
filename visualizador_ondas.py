#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  visualizador_ondas.py
#
#  Copyright 2025 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
import tkinter as tk
from tkinter import ttk
import numpy as np

# --- Importaciones de Matplotlib para Tkinter ---
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# --- MODELO (Lógica Matemática - SRP) ---
class WaveModel:
    """
    Se encarga únicamente de los cálculos matemáticos.
    No sabe nada de la interfaz gráfica.
    """
    def __init__(self):
        # Generamos el eje X (tiempo) una sola vez
        # De 0 a 4pi, con 500 puntos de resolución
        self.x = np.linspace(0, 4 * np.pi, 500)

    def calculate_wave(self, amp1, freq1, amp2, freq2):
        """
        Calcula la onda resultante combinada.
        Formula: y = A1*cos(f1*x) + A2*cos(f2*x)
        """
        y1 = amp1 * np.cos(freq1 * self.x)
        y2 = amp2 * np.cos(freq2 * self.x)
        y_total = y1 + y2
        return self.x, y_total, y1, y2

# --- VISTA / CONTROLADOR (Interfaz Gráfica - SRP) ---
class WaveApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Visualizador de Ondas Coseno")
        self.root.geometry("900x600")

        # Instancia del modelo
        self.model = WaveModel()

        # --- Configuración del Layout Principal ---
        # Dividimos la ventana en dos: Gráfico (Izquierda) y Controles (Derecha)
        self.frame_graph = tk.Frame(self.root)
        self.frame_graph.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.frame_controls = tk.Frame(self.root, bg="#f0f0f0", width=250)
        self.frame_controls.pack(side=tk.RIGHT, fill=tk.Y)

        # --- Inicialización ---
        self._init_graph()
        self._init_controls()

        # Dibujar por primera vez
        self.update_plot()

    def _init_graph(self):
        """Configura el lienzo de Matplotlib dentro de Tkinter."""
        # 1. Crear la Figura (el contenedor del gráfico)
        self.fig = Figure(figsize=(5, 4), dpi=100)

        # 2. Añadir un subplot (el área de dibujo)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title("Combinación de Ondas")
        self.ax.set_xlabel("Tiempo")
        self.ax.set_ylabel("Amplitud")
        self.ax.grid(True)

        # 3. Crear el Canvas especial de Tkinter
        # Este es el "puente" entre Matplotlib y Tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame_graph)
        self.canvas.draw()

        # 4. Empaquetar el widget del canvas en la interfaz
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def _init_controls(self):
        """Crea los deslizadores (Sliders) para controlar los parámetros."""

        lbl_title = tk.Label(self.frame_controls, text="Parámetros",
                             font=("Arial", 14, "bold"), bg="#f0f0f0")
        lbl_title.pack(pady=10)

        # --- Onda 1 ---
        self.slider_amp1 = self._create_slider("Amplitud 1", 0, 5, 2.0)
        self.slider_freq1 = self._create_slider("Frecuencia 1", 0.1, 5.0, 1.0)

        # Separador visual
        ttk.Separator(self.frame_controls, orient='horizontal').pack(fill='x', pady=15)

        # --- Onda 2 ---
        self.slider_amp2 = self._create_slider("Amplitud 2", 0, 5, 1.5)
        self.slider_freq2 = self._create_slider("Frecuencia 2", 0.1, 5.0, 3.0)

    def _create_slider(self, label_text, min_val, max_val, default):
        """Helper para crear sliders repetitivos (DRY - Don't Repeat Yourself)."""

        # Etiqueta
        lbl = tk.Label(self.frame_controls, text=label_text, bg="#f0f0f0")
        lbl.pack(pady=(10, 0))

        # Scale (Slider)
        # 'command' llama a self.update_plot cada vez que mueves el slider
        slider = tk.Scale(self.frame_controls, from_=min_val, to=max_val,
                          orient=tk.HORIZONTAL, resolution=0.1,
                          length=200, command=self.update_plot)
        slider.set(default)
        slider.pack()
        return slider

    def update_plot(self, event=None):
        """
        Callback: Se ejecuta cuando se mueve cualquier slider.
        Recalcula los datos y redibuja el gráfico.
        """
        # 1. Obtener valores actuales de los sliders
        a1 = self.slider_amp1.get()
        f1 = self.slider_freq1.get()
        a2 = self.slider_amp2.get()
        f2 = self.slider_freq2.get()

        # 2. Calcular datos usando el Modelo
        x, y_total, y1, y2 = self.model.calculate_wave(a1, f1, a2, f2)

        # 3. Limpiar y redibujar el gráfico
        # (Nota: Para máximo rendimiento en tiempo real se usa set_ydata,
        # pero clear() + plot() es más fácil de entender para empezar).
        self.ax.clear()

        # Dibujamos las ondas individuales (punteadas y finas)
        self.ax.plot(x, y1, color='green', alpha=0.3, linestyle='--', label='Onda 1')
        self.ax.plot(x, y2, color='blue', alpha=0.3, linestyle='--', label='Onda 2')

        # Dibujamos la suma (línea gruesa roja)
        self.ax.plot(x, y_total, color='red', linewidth=2, label='Suma')

        # Decoración del gráfico
        self.ax.legend(loc='upper right')
        self.ax.grid(True)
        self.ax.set_ylim(-10, 10) # Fijamos el eje Y para que no "salte"
        self.ax.set_title(f"Suma: {a1}cos({f1}x) + {a2}cos({f2}x)")

        # 4. Ordenar al canvas que se actualice
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = WaveApp(root)
    root.mainloop()
