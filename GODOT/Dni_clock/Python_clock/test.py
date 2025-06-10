import tkinter as tk
import datetime

def actualizar_hora():
    ahora = datetime.datetime.now()
    Year = ahora.year
    Month = ahora.month
    Day = ahora.day
    Hour = ahora.hour
    Minute = ahora.minute
    Second = ahora.second

    etiqueta_fecha.config(text=f"Fecha: {Year}-{Month:02}-{Day:02}")
    etiqueta_hora.config(text=f"Hora: {Hour:02}:{Minute:02}:{Second:02}")

    ventana.after(1000, actualizar_hora) # Actualizar cada segundo

ventana = tk.Tk()
ventana.title("Fecha y Hora")

# Etiquetas para la fecha y la hora
etiqueta_fecha = tk.Label(ventana, text="")
etiqueta_fecha.pack(pady=10)

etiqueta_hora = tk.Label(ventana, text="")
etiqueta_hora.pack(pady=10)

# Llamar a la función para actualizar la hora cada segundo
actualizar_hora()

ventana.mainloop()
