#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  File_organizer_tk.py
#
#  Copyright 2026 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
# ---------------------------
# File organizer, modificado para evitar sobrescribir y preguntar carpetas.
import os
import string
import shutil
#  1. Importar tkinter para los diálogos gráficos (GUI)
# Esta biblioteca viene incluida con Python, no necesita instalación.
from tkinter import Tk, filedialog

# --- Definiciones Globales (Constantes) ---

#  Eliminadas las rutas 'source_folder' y 'destination_root' hardcodeadas.

# Las definiciones de carpetas ahora son nombres relativos
DESTINATION_NAMES = {
    "Images": "Images",
    "Documents": "Documents",
    "Spreadsheets": "Spreadsheets",
    "Audio": "Audio",
    "Videos": "Videos",
    "Shortcuts": "Shortcuts",
    "Books": "Books",
    "Scripts": "Scripts",
    "Miscellaneous": "Miscellaneous",
    "Zipped": "Zipped",
    "Folders": "Folders",
}

# File type extensions for categories (Sin cambios)
FILE_CATEGORIES = {
    "Images": [".jpg", ".jpeg", ".png", ".gif", ".jfif", ".bmp", ".tiff", ".svg", ".html", ".webp"],
    "Documents": [".doc", ".docx", ".pdf", ".txt", ".ppt", ".pptx"],
    "Spreadsheets": [".xls", ".xlsx", ".csv", ".json", ".xml"],
    "Audio": [".mp3", ".wav", ".aac", ".flac", ".ogg", ".wma"],
    "Videos": [".mp4", ".mkv", ".avi", ".mov", ".wmv", ".flv"],
    "Scripts": [".ps1", ".py"],
    "Shortcuts": [".lnk", ".url", ".webloc"],
    "Zipped": [".zip", ".rar", ".7z"],
    "Books": [".epub"]
}

# --- Funciones de Lógica ---

def get_unique_path(destination_folder: str, item_name: str) -> str | None:
    """
    Genera una ruta de destino única para evitar sobrescribir.
    (Función sin cambios)
    """
    # 1. Dividir el nombre (funciona para "file.txt" y "Folder")
    base_name, extension = os.path.splitext(item_name)

    # 2. Comprobar la ruta original primero (KISS)
    target_path = os.path.join(destination_folder, item_name)
    if not os.path.exists(target_path):
        return target_path

    # 3. Si existe, iniciar la lógica de sufijos
    # Genera letras de 'a' a 'z'
    for char in string.ascii_lowercase:
        # Genera números de '01' a '99'
        for i in range(1, 100):
            num_str = f"{i:02d}"  # Formatea como '01', '02', ... '99'
            suffix = f"_({char}{num_str})"

            new_name = f"{base_name}{suffix}{extension}"
            new_target_path = os.path.join(destination_folder, new_name)

            if not os.path.exists(new_target_path):
                return new_target_path

    # Si se agotan todas las 2574 combinaciones (a01-z99)
    print(f"ERROR: No se pudo encontrar un nombre único para '{item_name}' en {destination_folder}")
    return None

def organize_files(source: str, destinations: dict):
    """
    Itera sobre la carpeta de origen, categoriza y mueve los items de forma segura.
     Modificado para aceptar 'destinations' como argumento (DIP)
    en lugar de usar una variable global.
    """
    for item in os.listdir(source):
        item_path = os.path.join(source, item)
        target_folder = None

        try:
            # 1. Determinar la carpeta de destino
            if os.path.isdir(item_path):
                #  Usa el diccionario 'destinations' inyectado
                target_folder = destinations["Folders"]

            elif os.path.isfile(item_path):
                _, extension = os.path.splitext(item)
                extension = extension.lower()

                # Buscar en categorías
                for category, extensions in FILE_CATEGORIES.items():
                    if extension in extensions:
                        #  Usa el diccionario 'destinations' inyectado
                        target_folder = destinations[category]
                        break
                else:
                    # Si no hay coincidencias, va a Misceláneos
                    #  Usa el diccionario 'destinations' inyectado
                    target_folder = destinations["Miscellaneous"]

            # 2. Si se encontró un destino, mover de forma segura
            if target_folder:
                # 2a. Obtener la ruta de destino ÚNICA
                unique_dest_path = get_unique_path(target_folder, item)

                # 2b. Mover el archivo si se encontró una ruta
                if unique_dest_path:
                    shutil.move(item_path, unique_dest_path)
                    print(f"Movido: {item} \n    -> {unique_dest_path}")
                else:
                    print(f"Omitido (Error de nombre): {item}")
            else:
                print(f"Omitido (Item no es archivo ni carpeta): {item}")

        except (IOError, OSError) as e:
            print(f"Error al mover {item}: {e}")
        except Exception as e:
            print(f"Error inesperado procesando {item}: {e}")

#  2. Nueva función para selección de carpetas
def get_folders_interactively() -> tuple[str | None, str | None]:
    """
    Usa un diálogo gráfico (Tkinter) para preguntar al usuario por las carpetas.

    Aplica SRP: Su única responsabilidad es la selección de carpetas.
    Retorna (ruta_origen, ruta_destino_raiz) o (None, None) si se cancela.
    """
    # Ocultar la ventana raíz de Tkinter
    root = Tk()
    root.withdraw()

    print("Iniciando selección de carpetas...")
    print("Por favor, selecciona la carpeta de ORIGEN (la que quieres organizar)...")
    source = filedialog.askdirectory(title="Selecciona la carpeta de ORIGEN")

    if not source:
        print("Operación cancelada. No se seleccionó carpeta de origen.")
        return None, None

    print(f"Carpeta de origen: {source}")
    print("Ahora, selecciona la carpeta de DESTINO (donde irá todo)...")
    destination_root = filedialog.askdirectory(title="Selecciona la carpeta de DESTINO RAÍZ")

    if not destination_root:
        print("Operación cancelada. No se seleccionó carpeta de destino.")
        return None, None

    print(f"Carpeta de destino raíz: {destination_root}")
    return source, destination_root

#  3. Función principal para encapsular la lógica
def main():
    """
    Función principal (wrapper) que ejecuta la aplicación.
    """
    # 1. Obtener carpetas del usuario
    source_folder, destination_root = get_folders_interactively()

    if not source_folder or not destination_root:
        print("Saliendo del programa.")
        return

    # 2. Construir el diccionario de destinos (DIP)
    #  Esto se mueve de global a aquí, después de la selección.
    print("Construyendo árbol de carpetas de destino...")
    destinations = {}
    for name, folder in DESTINATION_NAMES.items():
        destinations[name] = os.path.join(destination_root, folder)

    # 3. Crear carpetas de destino
    for folder_path in destinations.values():
        os.makedirs(folder_path, exist_ok=True)

    # 4. Ejecutar la organización
    print(f"\nIniciando organización de '{source_folder}'...")
    #  Inyectamos 'destinations' en la función
    organize_files(source_folder, destinations)
    print("\n¡Organización completada!")

if __name__ == "__main__":
    # El bloque __main__ ahora es mucho más limpio (KISS)
    main()
