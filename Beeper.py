import platform
import sys
import time

def hacer_beep(frecuencia=1000, duracion=500):
    """
    Emite un tono 'beep'.
    - frecuencia: En Hertzios (solo aplicable en Windows).
    - duracion: En milisegundos (solo aplicable en Windows).
    """
    sistema = platform.system()

    if sistema == "Windows":
        import winsound
        winsound.Beep(frecuencia, duracion)

    elif sistema == "Darwin":  # macOS
        # Utiliza la campana del sistema a través de la terminal
        sys.stdout.write('\a')
        sys.stdout.flush()

    else:  # Linux y otros Unix
        # Imprime el carácter ASCII Bell ('\a')
        sys.stdout.write('\a')
        sys.stdout.flush()
        # Espera para simular la duración si se llama consecutivamente
        time.sleep(duracion / 1000.0)

# Ejemplo de uso:
if __name__ == "__main__":
    print("Emitiendo beep...")
    hacer_beep(frecuencia=1000, duracion=500)
