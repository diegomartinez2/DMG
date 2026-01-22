#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
#
#  Copyright 2024 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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

# -------
# Clases
# -------
import cv2
from PIL import Image
import pytesseract
import speech_recognition as sr
import time
import sys
import os
import argparse

# PARÁMETRO CONFIGURABLE: Tiempo de grabación en segundos (default: 5)
TIEMPO_GARBADO = 5  # CAMBIA AQUÍ si quieres fijo (ej: 10, 3, 15)

# Configuración AUTOMÁTICA para Linux
pytesseract.pytesseract.tesseract_cmd = '/usr/bin/tesseract'

# Parser para argumentos de línea de comandos
parser = argparse.ArgumentParser(description='Captura imagen + audio con tiempo configurable')
parser.add_argument('--tiempo', '-t', type=int, default=TIEMPO_GARBADO,
                    help=f'Tiempo de grabación en segundos (default: {TIEMPO_GARBADO})')
args = parser.parse_args()

# Usa el tiempo del parámetro
tiempo_grabacion = args.tiempo
print(f"⏱️  Tiempo de grabación configurado: {tiempo_grabacion} segundos")

# Parte 1: Captura imagen
print("📷 Capturando imagen... (mira a la cámara)")
cap = cv2.VideoCapture(0)
ret, frame = cap.read()
if ret:
    cv2.imwrite('foto.png', frame)
cap.release()
print("✅ Imagen capturada")

descripcion_imagen = pytesseract.image_to_string(Image.open('foto.png'), lang='spa').strip()
if not descripcion_imagen:
    descripcion_imagen = "No se detectó texto en la imagen."

# Parte 2: Captura audio CON CONTADOR VISUAL DINÁMICO
print("\n🎤 Preparando micrófono...")
r = sr.Recognizer()
with sr.Microphone() as source:
    print("Habla ahora...")

    # CONTADOR VISUAL: Puntos que desaparecen según el tiempo
    puntos_totales = 5  # Siempre 5 puntos máximo, se ajustan proporcionalmente
    for i in range(tiempo_grabacion, 0, -1):
        puntos_restantes = int((i / tiempo_grabacion) * puntos_totales)
        sys.stdout.write(f"\r🎤 Grabando... {'·' * puntos_restantes} ({i}s)")
        sys.stdout.flush()
        time.sleep(1)

    print("\r🎤 Grabando...     ")  # Limpia la línea
    audio = r.listen(source, timeout=tiempo_grabacion)

try:
    transcripcion_audio = r.recognize_google(audio, language="es-ES")
except sr.UnknownValueError:
    transcripcion_audio = "No se entendió el audio."
except sr.RequestError:
    transcripcion_audio = "Error con el servicio."

# Formato de salida
output = f"[imagen] {descripcion_imagen} [imagen];[audio] {transcripcion_audio} [audio]"
print(f"\n📋 RESULTADO FINAL:\n{output}")

# Notificación nativa
os.system(f'notify-send "Captura Lista" "{output[:100]}..."')

# Limpieza
os.remove('foto.png')
print("🧹 Archivos temporales eliminados")
