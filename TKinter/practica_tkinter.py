#!/usr/bin/env python

import tkinter as tk
from tkinter import ttk, messagebox

# =====================================================================
# FUNCIONES DE CONTROL (Lógica simple)
# =====================================================================
def saludar():
    nombre = entrada_nombre.get()
    if nombre.strip():
        etiqueta_resultado.config(text=f"¡Hola, {nombre}! Bienvenido/a.")
    else:
        messagebox.showwarning("Atención", "Por favor, escribe un nombre.")

def actualizar_opciones():
    # Obtener el estado de los Checkbuttons
    opciones = []
    if var_check1.get(): opciones.append("Opción A")
    if var_check2.get(): opciones.append("Opción B")

    # Obtener el estado del Radiobutton
    seleccion_radio = var_radio.get()

    texto_resumen = f"Checks: {', '.join(opciones) if opciones else 'Ninguno'}\nRadio: {seleccion_radio}"
    etiqueta_resumen.config(text=texto_resumen)

# =====================================================================
# CONFIGURACIÓN DE LA VENTANA PRINCIPAL
# =====================================================================
ventana = tk.Tk()
ventana.title("Práctica de Tkinter - Filosofía KISS")
ventana.geometry("450x350")

# Creamos el contenedor de pestañas (Notebook)
panel_pestanas = ttk.Notebook(ventana)
panel_pestanas.pack(fill="both", expand=True, padx=10, pady=10)

# =====================================================================
# PESTAÑA 1: TEXTO Y ENTRADAS (Conceptos básicos)
# =====================================================================
pestana1 = ttk.Frame(panel_pestanas)
panel_pestanas.add(pestana1, text="Texto y Entradas")

# Componentes de la Pestaña 1
lbl_instruccion = tk.Label(pestana1, text="Escribe tu nombre abajo:", font=("Arial", 11))
lbl_instruccion.pack(pady=10)

entrada_nombre = tk.Entry(pestana1, font=("Arial", 11), width=25)
entrada_nombre.pack(pady=5)

btn_saludar = tk.Button(pestana1, text="Saludar", command=saludar, bg="#4CAF50", fg="white")
btn_saludar.pack(pady=10)

etiqueta_resultado = tk.Label(pestana1, text="", font=("Arial", 11, "bold"), fg="blue")
etiqueta_resultado.pack(pady=10)


# =====================================================================
# PESTAÑA 2: SELECCIONES (Checks y Radios)
# =====================================================================
pestana2 = ttk.Frame(panel_pestanas)
panel_pestanas.add(pestana2, text="Selecciones")

# Variables de control para almacenar los estados
var_check1 = tk.BooleanVar()
var_check2 = tk.BooleanVar()
var_radio = tk.StringVar(value="No seleccionado")

# Grupo de Checkbuttons (Selección múltiple)
lbl_check = tk.Label(pestana2, text="Selección Múltiple (Checkbuttons):", font=("Arial", 10, "bold"))
lbl_check.pack(anchor="w", padx=20, pady=5)

ch_opcion1 = tk.Checkbutton(pestana2, text="Activar Opción A", variable=var_check1, command=actualizar_opciones)
ch_opcion1.pack(anchor="w", padx=40)

ch_opcion2 = tk.Checkbutton(pestana2, text="Activar Opción B", variable=var_check2, command=actualizar_opciones)
ch_opcion2.pack(anchor="w", padx=40)

# Separador visual simple
tk.Label(pestana2, text="-"*40, fg="gray").pack(pady=5)

# Grupo de Radiobuttons (Selección única)
lbl_radio = tk.Label(pestana2, text="Selección Única (Radiobuttons):", font=("Arial", 10, "bold"))
lbl_radio.pack(anchor="w", padx=20, pady=5)

r_si = tk.Radiobutton(pestana2, text="Modo Claro", value="Claro", variable=var_radio, command=actualizar_opciones)
r_si.pack(anchor="w", padx=40)

r_no = tk.Radiobutton(pestana2, text="Modo Oscuro", value="Oscuro", variable=var_radio, command=actualizar_opciones)
r_no.pack(anchor="w", padx=40)

# Cuadro inferior de resumen en tiempo real
etiqueta_resumen = tk.Label(pestana2, text="Checks: Ninguno\nRadio: No seleccionado", fg="green", justify="left")
etiqueta_resumen.pack(pady=15)


# =====================================================================
# PESTAÑA 3: ALERTAS (Interacción con el sistema)
# =====================================================================
pestana3 = ttk.Frame(panel_pestanas)
panel_pestanas.add(pestana3, text="Alertas")

lbl_alertas = tk.Label(pestana3, text="Prueba los diferentes diálogos del sistema:", font=("Arial", 10))
lbl_alertas.pack(pady=15)

# Funciones inline rápidas (lambdas) para no sobrecargar el código de arriba
btn_info = tk.Button(pestana3, text="Mostrar Info",
                     command=lambda: messagebox.showinfo("Información", "Esto es un mensaje informativo."))
btn_info.pack(pady=5, fill="x", padx=50)

btn_error = tk.Button(pestana3, text="Mostrar Error",
                      command=lambda: messagebox.showerror("Error", "Algo salió mal (ejemplo)."))
btn_error.pack(pady=5, fill="x", padx=50)

btn_salir = tk.Button(pestana3, text="Cerrar Aplicación", bg="#f44336", fg="white",
                      command=lambda: ventana.quit() if messagebox.askyesno("Salir", "¿Seguro que quieres salir?") else None)
btn_salir.pack(pady=20, fill="x", padx=50)


# =====================================================================
# BUCLE PRINCIPAL
# =====================================================================
ventana.mainloop()
