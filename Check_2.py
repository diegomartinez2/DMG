#!/usr/bin/env python
import psutil
import os

def procesos_en_directorio(directorio):
    """
    Verifica procesos cuyo cwd es 'directorio'.
    Retorna: lista de (PID, nombre) o [] si ninguno.
    """
    dir_abs = os.path.abspath(directorio)
    en_directorio = []
    for proc in psutil.process_iter(['pid', 'name', 'cwd']):
        try:
            if proc.info['cwd'] == dir_abs:
                en_directorio.append((proc.info['pid'], proc.info['name']))
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return en_directorio

# Uso
resultado = procesos_en_directorio('.')
if resultado:
    print("Procesos corriendo en el directorio:")
    for pid, nombre in resultado:
        print(f"  PID {pid}: {nombre}")
else:
    print("No hay procesos con cwd en el directorio actual.")
