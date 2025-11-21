#!/usr/bin/env python
import tkinter as tk
from tkinter import messagebox
'''
Este es un ejemplo ilustrativo de cómo usar tkinter en un script python.
Aquí se muestran los conceptos básicos para su uso, los tres pilares fundamentales de tkinter:

    Widgets: Los bloques de construcción (botones, etiquetas, entradas).

    Geometría (Layout): Cómo organizar las cosas (usaremos .pack() y .grid()).

    Eventos: Cómo hacer que los botones hagan algo (callbacks).

Conceptos Clave Explicados

    tk.Tk() y .mainloop():

        root = tk.Tk() crea la ventana principal. Es el lienzo.

        root.mainloop() es un bucle infinito que espera clics o teclas. Sin esto, la ventana aparecería y se cerraría en un microsegundo.

    Los Managers de Geometría (pack vs grid):

        .pack(): Es el más simple (KISS). "Pon esto arriba, luego esto abajo". Útil para listas verticales simples.

        .grid(): Es como una hoja de cálculo (fila 0, columna 1). Te da mucho más control. En el ejemplo, lo usé para poner los botones de porcentaje uno al lado del otro.

    Widgets:

        Label: Muestra texto.

        Entry: Pide texto.

        Button: Ejecuta una acción.

        Frame: Agrupa otros widgets.

    Eventos (command=):

        Los botones tienen un argumento command. Le pasas el nombre de una función.

        Truco Pro: Si necesitas pasar argumentos (como el 10, 15 o 20 en mi ejemplo), usas lambda: funcion(10).
'''
# --- Lógica del Negocio (Separada de la UI - SRP) ---
class TipCalculatorModel:
    """
    Maneja únicamente la lógica matemática.
    No sabe nada sobre ventanas o botones.
    """
    def calculate_total(self, bill_amount: float, tip_percent: int) -> float:
        if bill_amount < 0:
            raise ValueError("La cuenta no puede ser negativa")

        tip_amount = bill_amount * (tip_percent / 100)
        return bill_amount + tip_amount

# --- Interfaz Gráfica (UI - SRP) ---
class TipCalculatorApp:
    """
    Maneja únicamente la interfaz visual y la interacción del usuario.
    """
    def __init__(self, root):
        self.root = root
        self.model = TipCalculatorModel() # Conexión con la lógica

        # 1. Configuración de la Ventana Principal
        self.root.title("Calculadora de Propinas")
        self.root.geometry("300x250") # Ancho x Alto

        # 2. Creación de Widgets (Elementos visuales)
        self._create_widgets()

    def _create_widgets(self):
        # --- Etiqueta y Entrada para el Total de la Cuenta ---
        # Label: Texto estático
        lbl_bill = tk.Label(self.root, text="Total de la Cuenta ($):")
        lbl_bill.pack(pady=5) # pack() coloca el elemento uno debajo del otro
                              # pady=5 añade 5 píxeles de espacio vertical

        # Entry: Campo de texto para escribir
        self.entry_bill = tk.Entry(self.root)
        self.entry_bill.pack(pady=5)

        # --- Botones para el Porcentaje (Usando Frame para agrupar) ---
        # Frame: Un contenedor invisible para organizar otros widgets
        frame_percentages = tk.Frame(self.root)
        frame_percentages.pack(pady=10)

        # Creamos 3 botones dentro del Frame
        # command=... le dice a Python qué función ejecutar al hacer clic
        # Usamos lambda para pasar argumentos a la función
        btn_10 = tk.Button(frame_percentages, text="10%",
                           command=lambda: self._on_calculate(10))
        btn_10.grid(row=0, column=0, padx=5) # grid() organiza en filas/columnas

        btn_15 = tk.Button(frame_percentages, text="15%",
                           command=lambda: self._on_calculate(15))
        btn_15.grid(row=0, column=1, padx=5)

        btn_20 = tk.Button(frame_percentages, text="20%",
                           command=lambda: self._on_calculate(20))
        btn_20.grid(row=0, column=2, padx=5)

        # --- Etiqueta para el Resultado ---
        self.lbl_result = tk.Label(self.root, text="Total a Pagar: ---",
                                   font=("Arial", 12, "bold"), fg="blue")
        self.lbl_result.pack(pady=20)

    def _on_calculate(self, percent):
        """
        Manejador del evento de clic (El puente entre UI y Lógica).
        """
        try:
            # 1. Obtener texto del input
            bill_text = self.entry_bill.get()

            # 2. Validar y convertir
            if not bill_text:
                return # No hacer nada si está vacío

            bill_amount = float(bill_text)

            # 3. Llamar a la lógica (Model)
            total = self.model.calculate_total(bill_amount, percent)

            # 4. Actualizar la UI
            self.lbl_result.config(text=f"Total a Pagar: ${total:.2f}")

        except ValueError:
            # Mostrar una ventana emergente de error nativa
            messagebox.showerror("Error", "Por favor, ingresa un número válido.")

# --- Punto de Entrada ---
if __name__ == "__main__":
    # 1. Crear la ventana raíz (root)
    root = tk.Tk()

    # 2. Instanciar nuestra App
    app = TipCalculatorApp(root)

    # 3. Iniciar el bucle principal (mantiene la ventana abierta)
    root.mainloop()
