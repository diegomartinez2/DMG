#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Check_1.py
#
#  Copyright 2022 Diego <diego@u038025>
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
