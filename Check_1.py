#!/usr/bin/env python
import psutil
import os

def archivos_en_uso(directorio):
    """
    Verifica si algún archivo en 'directorio' está abierto por un proceso.
    Retorna: lista de (PID, archivo) o [] si no hay nada en uso.
    """
    dir_abs = os.path.abspath(directorio)
    en_uso = []
    for proc in psutil.process_iter(['pid', 'name']):
        try:
            for file in proc.open_files():
                if file.path.startswith(dir_abs):
                    en_uso.append((proc.info['pid'], proc.info['name'], file.path))
        except (psutil.NoSuchProcess, psutil.AccessDenied, OSError):
            pass  # Ignora procesos inaccesibles
    return en_uso

# Uso
resultado = archivos_en_uso('.')
if resultado:
    print("Archivos en uso:")
    for pid, nombre, archivo in resultado:
        print(f"  PID {pid} ({nombre}): {archivo}")
else:
    print("No hay archivos abiertos en el directorio actual.")
