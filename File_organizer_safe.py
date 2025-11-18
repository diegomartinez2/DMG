#!/usr/bin/env python
# File organizer, modificado para evitar sobrescribir archivos.
import os
import string
import shutil

# Define source folder to analyze
# !!! ACTUALIZA ESTA RUTA a tu carpeta de origen
source_folder = r"C:\Users\USERNAME\Desktop\to_organize"

# Define destination folder paths
# !!! ACTUALIZA ESTA RUTA a tu carpeta de destino
destination_root = r"C:\Users\USERNAME\Desktop\Organized"
DESTINATIONS = {
    "Images": os.path.join(destination_root, "Images"),
    "Documents": os.path.join(destination_root, "Documents"),
    "Spreadsheets": os.path.join(destination_root, "Spreadsheets"),
    "Audio": os.path.join(destination_root, "Audio"),
    "Videos": os.path.join(destination_root, "Videos"),
    "Shortcuts": os.path.join(destination_root, "Shortcuts"),
    "Books": os.path.join(destination_root, "Books"),
    "Scripts": os.path.join(destination_root, "Scripts"),
    "Miscellaneous": os.path.join(destination_root, "Miscellaneous"),
    "Zipped": os.path.join(destination_root, "Zipped"),
    "Folders": os.path.join(destination_root, "Folders"),
}

# File type extensions for categories
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

# Create destination folders if they don't exist
for folder in DESTINATIONS.values():
    os.makedirs(folder, exist_ok=True)

def get_unique_path(destination_folder: str, item_name: str) -> str | None:
    """
    Genera una ruta de destino única para evitar sobrescribir.

    Aplica SRP: Esta función tiene la única responsabilidad de
    resolver colisiones de nombres.

    Sigue la lógica: "file.txt" -> "file_(a01).txt" -> ... -> "file_(z99).txt"

    Args:
        destination_folder (str): La carpeta donde se moverá el item.
        item_name (str): El nombre del archivo o carpeta.

    Returns:
        str | None: La ruta de destino única y segura, o None si falla.
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

def organize_files(source: str):
    """
    Itera sobre la carpeta de origen, categoriza y mueve los items de forma segura.
    """
    for item in os.listdir(source):
        item_path = os.path.join(source, item)
        target_folder = None

        try:
            # 1. Determinar la carpeta de destino
            if os.path.isdir(item_path):
                target_folder = DESTINATIONS["Folders"]

            elif os.path.isfile(item_path):
                _, extension = os.path.splitext(item)
                extension = extension.lower()

                # Buscar en categorías
                for category, extensions in FILE_CATEGORIES.items():
                    if extension in extensions:
                        target_folder = DESTINATIONS[category]
                        break
                else:
                    # Si no hay coincidencias, va a Misceláneos
                    target_folder = DESTINATIONS["Miscellaneous"]

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


if __name__ == "__main__":
    if not os.path.exists(source_folder):
        print(f"La carpeta de origen '{source_folder}' no existe. Por favor, revisa la ruta en la variable 'source_folder'.")
    elif not os.path.exists(destination_root):
        print(f"La carpeta de destino '{destination_root}' no existe. Se creará.")
        os.makedirs(destination_root, exist_ok=True)
        organize_files(source_folder)
        print("\nOrganización completada!")
    else:
        print("Iniciando organización...")
        organize_files(source_folder)
        print("\nOrganización completada!")
