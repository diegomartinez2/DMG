#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   morse_code_player.py
#
#   Copyright 2026 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#   MA 02110-1301, USA.
#

import curses
import time
import sys
import math
import struct
import subprocess
import threading

# --- Módulo 1: Traductor (SRP) ---
# Esta clase solo tiene una responsabilidad: traducir texto a Morse.

class MorseTranslator:
    """
    Traduce cadenas de texto de/hacia el código Morse Internacional.
    Aplica SRP: su única responsabilidad es la traducción.
    """

    MORSE_CODE_DICT = {
        'A': '.-', 'B': '-...', 'C': '-.-.', 'D': '-..', 'E': '.',
        'F': '..-.', 'G': '--.', 'H': '....', 'I': '..', 'J': '.---',
        'K': '-.-', 'L': '.-..', 'M': '--', 'N': '-.', 'O': '---',
        'P': '.--.', 'Q': '--.-', 'R': '.-.', 'S': '...', 'T': '-',
        'U': '..-', 'V': '...-', 'W': '.--', 'X': '-..-', 'Y': '-.--',
        'Z': '--..',
        '1': '.----', '2': '..---', '3': '...--', '4': '....-', '5': '.....',
        '6': '-....', '7': '--...', '8': '---..', '9': '----.', '0': '-----',
        ',': '--..--', '.': '.-.-.-', '?': '..--..', '/': '-..-.',
        '-': '-....-', '(': '-.--.', ')': '-.--.-',
        ' ': '/' # Usamos '/' para el espacio entre palabras
    }

    def get_static_display(self, text: str) -> str:
        """
        Traduce texto a una cadena estática de Morse para visualización.
        """
        morse_text = ""
        for char in text.upper():
            morse_char = self.MORSE_CODE_DICT.get(char, ' ')

            if morse_char == '/':
                morse_text += " /  "
            elif morse_char != ' ':
                morse_text += morse_char + " "

        return morse_text.strip()


# --- Módulo 2: Reproductor Visual y Sonoro (SRP/OCP) ---
# Esta clase maneja la E/S de la terminal (curses), la temporización
# y el audio vía aplay.

class MorsePlayer:
    """
    Maneja la reproducción visual y sonora del código Morse.
    Aplica SRP: se encarga de la visualización, temporización y emisión de tonos.
    """

    def __init__(self, stdscr, dot_time: float, dash_time: float, frequency: int = 800):
        """
        Inicializa el reproductor, la interfaz de curses y la frecuencia del tono.
        """
        self.stdscr = stdscr
        self.frequency = frequency

        # Configuración de Tiempos (OCP)
        T = dot_time
        self.T_DOT = T
        self.T_DASH = dash_time
        self.T_SYMBOL_SPACE = T      # 1T entre símbolos
        self.T_LETTER_SPACE = T * 3  # 3T entre letras
        self.T_WORD_SPACE = T * 7    # 7T entre palabras

        # Configuración de Curses
        curses.curs_set(0)
        self.stdscr.nodelay(True)
        self.stdscr.timeout(100)

        max_y, max_x = self.stdscr.getmaxyx()

        # Ventanas
        h_static = max_y // 2
        self.win_static = stdscr.subwin(h_static, max_x, 0, 0)
        self.win_static.box()
        self.win_static.addstr(1, 2, "Morse Estático (Puntos y Rayas):")

        h_blinker = max_y - h_static
        self.win_blinker = stdscr.subwin(h_blinker, max_x, h_static, 0)
        self.win_blinker.box()
        self.win_blinker.addstr(1, 2, "Transmisión (Asterisco y Sonido): [Presiona 'q' para salir]")

        self.stdscr.refresh()

    def _emitir_beep_async(self, duracion_seg: float):
        """
        Genera y reproduce el tono en segundo plano mediante un hilo,
        evitando bloquear el bucle principal de curses.
        """
        def _play():
            duracion_ms = duracion_seg * 1000.0
            sample_rate = 44100
            n_samples = int(sample_rate * (duracion_ms / 1000.0))

            audio_data = bytearray()
            for i in range(n_samples):
                val = int(32767 * 0.3 * math.sin(2 * math.pi * self.frequency * i / sample_rate))
                audio_data.extend(struct.pack('<h', val))

            try:
                proceso = subprocess.Popen(
                    ['aplay', '-q', '-f', 'S16_LE', '-r', str(sample_rate), '-c', '1'],
                    stdin=subprocess.PIPE
                )
                proceso.communicate(audio_data)
            except Exception:
                pass  # Ignora errores si aplay no está presente en el sistema

        threading.Thread(target=_play, daemon=True).start()

    def _update_blinker(self, text: str):
        """Actualiza la ventana de parpadeo visual."""
        self.win_blinker.clear()
        self.win_blinker.box()
        self.win_blinker.addstr(1, 2, "Transmisión (Asterisco y Sonido): [Presiona 'q' para salir]")

        h, w = self.win_blinker.getmaxyx()
        self.win_blinker.addstr(h // 2, w // 2 - len(text) // 2, text, curses.A_BOLD)
        self.win_blinker.refresh()

    def _sleep_with_exit(self, duration: float) -> bool:
        """Pausa la ejecución mientras comprueba si se presiona 'q'."""
        steps = int(duration / 0.05)
        if steps == 0:
            steps = 1

        for _ in range(steps):
            if self.stdscr.getch() == ord('q'):
                return True
            time.sleep(duration / steps)
        return False

    def play_text(self, text: str, translator: MorseTranslator):
        """Punto de entrada principal para reproducir un texto visual y sonoramente."""
        static_morse = translator.get_static_display(text)
        max_w = self.win_static.getmaxyx()[1] - 4
        self.win_static.addstr(3, 4, static_morse[:max_w])
        if len(static_morse) > max_w:
            self.win_static.addstr(4, 4, static_morse[max_w:max_w*2])
        self.win_static.refresh()

        for i, char in enumerate(text.upper()):
            morse_char = translator.MORSE_CODE_DICT.get(char)

            if morse_char == '/':
                self._update_blinker("  /  ")
                if self._sleep_with_exit(self.T_WORD_SPACE - self.T_LETTER_SPACE):
                    return

            elif morse_char:
                for j, symbol in enumerate(morse_char):
                    # A. Determinar duración y emitir SÍMBOLO (*) + AUDIO
                    duration = self.T_DOT if symbol == '.' else self.T_DASH
                    self._update_blinker("*")
                    self._emitir_beep_async(duration)

                    if self._sleep_with_exit(duration):
                        return

                    # B. Mostrar el ESPACIO entre símbolos (1T)
                    self._update_blinker(" ")
                    if self._sleep_with_exit(self.T_SYMBOL_SPACE):
                        return

                # C. Mostrar el ESPACIO entre letras (3T)
                if i < len(text) - 1 and text[i+1] != ' ':
                    if self._sleep_with_exit(self.T_LETTER_SPACE - self.T_SYMBOL_SPACE):
                        return

        self._update_blinker("[Transmisión completada. 'q' para salir.]")
        while self.stdscr.getch() != ord('q'):
            time.sleep(0.1)


# --- Módulo 3: Aplicación Principal (OCP/DIP) ---
'''
def main_app(stdscr, text_to_play: str):
    DOT_DURATION = 0.25   # Duración del punto (0.25s)
    DASH_DURATION = 1.0   # Duración de la raya (1.0s)
    FREQ_HZ = 800         # Frecuencia del sonido en Hz
'''
def main_app(stdscr, text_to_play: str):
    # Definir velocidad en WPM directamente
    WPM = 15  # Cambia este valor (ej. 12, 15, 20)

    DOT_DURATION = 1.2 / WPM
    DASH_DURATION = DOT_DURATION * 3

    FREQ_HZ = 800

    translator = MorseTranslator()
    player = MorsePlayer(stdscr, DOT_DURATION, DASH_DURATION, FREQ_HZ)

    try:
        player.play_text(text_to_play, translator)
    except curses.error:
        pass
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    texto_ejemplo = "SOS HABI TEST"

    if len(sys.argv) > 1:
        texto_ejemplo = " ".join(sys.argv[1:])

    print(f"Iniciando transmisor Morse para: '{texto_ejemplo}'")
    print("La aplicación tomará control de la terminal...")
    print("Presiona 'q' en cualquier momento para salir.")
    time.sleep(2.5)

    curses.wrapper(main_app, texto_ejemplo)

    print("Transmisor detenido. Terminal restaurada.")
