# importar la libreria TK
import tkinter as tk
from tkinter import ttk
# nota las funciones han de situarse antes del código que las llama
def function_click():
    accion.configure(text="Has pusado el boton")
    etiqueta.configure(foreground='red')
# iniciar una ventana
ventana = tk.Tk()
ventana.title("python-Tkinter")
# para evitar el redimensionado de la ventana
#ventana.resizable(0,0)
# Agregar una etiqueta
ttk.Label(ventana, text="Texto de la etiqueta").grid(column=0,row=0)
# Otra forma de agregar etiquetas
etiqueta=ttk.Label(ventana,text="Otra etiqueta")
etiqueta.grid(column=0,row=1)

# Agregar un Boton
accion=ttk.Button(ventana, text="pulsa", command=function_click)
accion.grid(column=1,row=0)



# Activar la ventana
ventana.mainloop()
