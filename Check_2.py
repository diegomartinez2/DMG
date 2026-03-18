#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Check_2.py
#
#  Copyright 2026 Diego <diego@u038025>
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
