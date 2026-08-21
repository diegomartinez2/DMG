import math
import struct
import subprocess
import time

def hacer_beep(duracion_ms, frecuencia=800):
    """Genera un tono puro y lo reproduce en Linux a través de aplay."""
    sample_rate = 44100
    n_samples = int(sample_rate * (duracion_ms / 1000.0))

    # Generar la onda senoidal en PCM de 16 bits
    audio_data = bytearray()
    for i in range(n_samples):
        # Escalar amplitud a un nivel cómodo (0.3 de volumen máximo)
        val = int(32767 * 0.3 * math.sin(2 * math.pi * frecuencia * i / sample_rate))
        audio_data.extend(struct.pack('<h', val))

    # Enviar los datos raw directamente a aplay
    proceso = subprocess.Popen(
        ['aplay', '-q', '-f', 'S16_LE', '-r', str(sample_rate), '-c', '1'],
        stdin=subprocess.PIPE
    )
    proceso.communicate(audio_data)

# Funciones de Código Morse
def punto():
    hacer_beep(duracion_ms=100)  # Corto
    time.sleep(0.1)             # Pausa entre elementos

def raya():
    hacer_beep(duracion_ms=300)  # Largo (3 veces el punto)
    time.sleep(0.1)

def pausa_letra():
    time.sleep(0.3)

# Ejemplo: Enviar la señal SOS (... --- ...)
if __name__ == "__main__":
    print("Enviando S.O.S. en código Morse...")

    # S (...)
    punto(); punto(); punto()
    pausa_letra()

    # O (---)
    raya(); raya(); raya()
    pausa_letra()

    # S (...)
    punto(); punto(); punto()
