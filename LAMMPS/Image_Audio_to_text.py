import cv2
from PIL import Image
import pytesseract
import speech_recognition as sr
import time
import sys
import os

# Configuración AUTOMÁTICA para Linux (Ubuntu/Debian)
pytesseract.pytesseract.tesseract_cmd = '/usr/bin/tesseract'  # Path estándar en Linux

# Parte 1: Captura imagen de cámara y extrae texto
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

# Parte 2: Captura audio CON CONTADOR VISUAL
print("\n🎤 Preparando micrófono...")
r = sr.Recognizer()
with sr.Microphone() as source:
    print("Habla ahora...")

    # CONTADOR VISUAL: 5 puntos que desaparecen
    for i in range(5, 0, -1):
        sys.stdout.write(f"\r🎤 Grabando... {'.' * i} ({i}s)")
        sys.stdout.flush()
        time.sleep(1)

    print("\r🎤 Grabando...     ")  # Limpia la línea
    audio = r.listen(source, timeout=5)

try:
    transcripcion_audio = r.recognize_google(audio, language="es-ES")
except sr.UnknownValueError:
    transcripcion_audio = "No se entendió el audio."
except sr.RequestError:
    transcripcion_audio = "Error con el servicio."

# Formato de salida
output = f"[imagen] {descripcion_imagen} [imagen];[audio] {transcripcion_audio} [audio]"
print(f"\n📋 RESULTADO FINAL:\n{output}")

# Notificación nativa de Linux (¡aparece en el escritorio!)
os.system(f'notify-send "Captura Lista" "{output[:100]}..."')

# Limpieza
os.remove('foto.png')
print("🧹 Archivos temporales eliminados")
