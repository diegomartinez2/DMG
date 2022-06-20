# importar la libreria TK
import tkinter as tk
from tkinter import ttk
# iniciar una ventana
ventana = tk.Tk()
ventana.title("python-Tkinter")
# para evitar el redimensionado de la ventana
#ventana.resizable(0,0)
# Agregar una etiqueta
ttk.Label(ventana, text="Texto de la etiqueta").grid(column=0,row=0)

# Activar la ventana
ventana.mainloop()
