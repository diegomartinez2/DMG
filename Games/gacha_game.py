#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  gacha_game.py
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
# ---------------------------
# Importación de los módulos
# ---------------------------import random
import tkinter as tk
from tkinter import ttk

# Función para cargar objetos desde un archivo de texto
def cargar_objetos(archivo):
    objetos = {}
    try:
        with open(archivo, 'r') as f:
            for linea in f:
                nombre, probabilidad = linea.strip().split(',')
                try:
                    probabilidad = float(probabilidad)
                    objetos[nombre] = Objeto(nombre, probabilidad)
                except ValueError:
                    print(f"Error al convertir la probabilidad '{probabilidad}' para el objeto '{nombre}'. Ignorando esta línea.")
    except FileNotFoundError:
        print(f"Archivo '{archivo}' no encontrado.  Se utilizarán objetos por defecto.")
        #  Si el archivo no existe, se utiliza un conjunto de objetos por defecto
        objetos = {
            Objeto("Espada de Hierro", 0.3),
            Objeto("Escudo de Madera", 0.2),
            Objeto("Poción de Curación", 0.4),
            Objeto("Amuleto de la Buena Fortuna", 0.1)
        }
    return objetos


# Clase para representar un objeto del gacha
class Objeto:
    def __init__(self, nombre, probabilidad):
        self.nombre = nombre
        self.probabilidad = probabilidad

    def __repr__(self):
        return f"Objeto: {self.nombre} (Probabilidad: {self.probabilidad})"

# Clase para el Gacha
class Gacha:
    def __init__(self, objetos):
        self.objetos = objetos  # Diccionario de objetos

    def sacar(self):
        # Implementación basada en probabilidades
        objeto_sacado = random.choices(list(self.objetos.keys()), weights=list(self.objetos.values()))[0]
        return objeto_sacado

# Clase para la interfaz de usuario (GUI)
class GachaGUI:
    def __init__(self, gacha):
        self.gacha = gacha
        self.window = tk.Tk()
        self.window.title("Juego Gacha")

        # Crear un frame para el inventario
        self.inventory_frame = ttk.Frame(self.window)
        self.inventory_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Etiqueta para el inventario
        self.inventory_label = ttk.Label(self.inventory_frame, text="Inventario:")
        self.inventory_label.pack(pady=5)

        # Crear un listbox para mostrar el inventario
        self.inventory_listbox = tk.Listbox(self.inventory_frame, width=30, height=10)
        self.inventory_listbox.pack(fill=tk.BOTH, expand=True)
        self.actualizar_inventario()

        # Botón para sacar un objeto
        self.sacar_button = tk.Button(self.window, text="Sacar Objeto", command=self.actualizar_inventario)
        self.sacar_button.pack(pady=10)

    def actualizar_inventario(self):
        # Simula la extracción de un objeto
        objeto_sacado = self.gacha.sacar()
        self.inventory_listbox.insert(tk.END, str(objeto_sacado))
        self.inventory_listbox.see(tk.END) # Scroll al final del listbox

    def run(self):
        self.window.mainloop()


# Configuración inicial
if __name__ == "__main__":
    archivo_objetos = "objetos_gacha.txt"  # Nombre del archivo de texto
    # Ejemplo de contenido del archivo "objetos_gacha.txt":
    # Espada de Hierro,0.3
    # Escudo de Madera,0.2
    # Poción de Curación,0.4
    # Amuleto de la Buena Fortuna,0.1

    objetos_gacha = cargar_objetos(archivo_objetos)

    # Crear instancia del Gacha
    gacha = Gacha(objetos_gacha)

    # Crear instancia de la GUI
    gui = GachaGUI(gacha)

    # Ejecutar la GUI
    gui.run()
