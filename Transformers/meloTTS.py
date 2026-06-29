'''
sudo apt-get update && sudo apt-get install -y python3-pip python3-venv espeak-ng ffmpeg git

git clone https://github.com/myshell-ai/MeloTTS.git
cd MeloTTS
pip install -e .

sudo apt-get install -y portaudio19-dev
pip install sounddevice numpy
'''

import os
import numpy as np
import sounddevice as sd
from melo.api import TTS

# 1. Configurar MeloTTS
model = TTS(language='ES', device='cpu') # Cambia a 'cuda:0' si usas GPU
speaker_id = model.hps.data.spk2id['ES']
sampling_rate = model.hps.data.sampling_rate

# 2. Leer el contenido del archivo de texto
archivo_texto = "texto.txt"

if not os.path.exists(archivo_texto):
    print(f"Error: El archivo '{archivo_texto}' no existe. Créalo primero.")
    exit()

with open(archivo_texto, 'r', encoding='utf-8') as f:
    texto_completo = f.read()

print("Generando audio y reproduciendo...")

# 3. Generar el audio directamente en memoria (retorna un array de NumPy)
audio_array = model.tts_to_file(texto_completo, speaker_id, speed=1.0)

# 4. Reproducir directamente por los altavoces
sd.play(audio_array, sampling_rate)
sd.wait() # Espera a que termine la reproducción antes de cerrar el script

print("Reproducción finalizada.")
