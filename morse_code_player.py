import curses
import time
import sys

# --- Módulo 1: Traductor (SRP) ---
# Esta clase solo tiene una responsabilidad: traducir texto a Morse.
# Es reutilizable en cualquier otro script (web, GUI, etc.)

class MorseTranslator:
    """
    Traduce cadenas de texto de/hacia el código Morse Internacional.
    Aplica SRP: su única responsabilidad es la traducción.
    """

    # El diccionario de traducción
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

        Args:
            text (str): El texto a traducir.

        Returns:
            str: Una cadena formateada con código Morse.
        """
        morse_text = ""
        for char in text.upper():
            morse_char = self.MORSE_CODE_DICT.get(char, ' ') # Obtener morse o un espacio

            if morse_char == '/':
                morse_text += " /  " # Separador de palabras
            elif morse_char != ' ':
                morse_text += morse_char + " " # Separador de letras

        return morse_text.strip()


# --- Módulo 2: Reproductor (SRP) ---
# Esta clase maneja toda la E/S de la terminal (curses), el diseño
# y la temporización. No sabe nada de la lógica de traducción.

class MorsePlayer:
    """
    Maneja la reproducción visual del código Morse usando curses.

    Aplica SRP: su única responsabilidad es la visualización y temporización.
    Aplica OCP: los tiempos son configurables en el constructor.
    """

    def __init__(self, stdscr, dot_time: float, dash_time: float):
        """
        Inicializa el reproductor y la interfaz de curses.

        Args:
            stdscr: La ventana principal de curses.
            dot_time (float): Duración (segundos) para un punto.
            dash_time (float): Duración (segundos) para una raya.
        """
        self.stdscr = stdscr

        # --- Configuración de Tiempos (OCP) ---
        # Tiempos base (1T = dot_time)
        T = dot_time
        self.T_DOT = T
        # El estándar es 3T, pero usamos el solicitado
        self.T_DASH = dash_time
        self.T_SYMBOL_SPACE = T      # 1T de espacio entre símbolos
        self.T_LETTER_SPACE = T * 3  # 3T de espacio entre letras
        self.T_WORD_SPACE = T * 7    # 7T de espacio entre palabras

        # --- Configuración de Curses ---
        curses.curs_set(0)  # Ocultar el cursor
        self.stdscr.nodelay(True) # No bloquear en getch()
        self.stdscr.timeout(100)  # Revisar entrada cada 100ms

        # Obtener dimensiones
        max_y, max_x = self.stdscr.getmaxyx()

        # --- Creación de ventanas (Pantalla dividida) ---
        h_static = max_y // 2
        self.win_static = stdscr.subwin(h_static, max_x, 0, 0)
        self.win_static.box()
        self.win_static.addstr(1, 2, "Morse Estático (Puntos y Rayas):")

        h_blinker = max_y - h_static
        self.win_blinker = stdscr.subwin(h_blinker, max_x, h_static, 0)
        self.win_blinker.box()
        self.win_blinker.addstr(1, 2, "Transmisión (Asterisco): [Presiona 'q' para salir]")

        self.stdscr.refresh()

    def _update_blinker(self, text: str):
        """
        Función interna (SRP) para actualizar solo la ventana de parpadeo.
        (KISS: Mantiene la lógica de dibujo simple y en un solo lugar).
        """
        self.win_blinker.clear()
        self.win_blinker.box()
        self.win_blinker.addstr(1, 2, "Transmisión (Asterisco): [Presiona 'q' para salir]")

        h, w = self.win_blinker.getmaxyx()
        # Centrar el texto
        self.win_blinker.addstr(h // 2, w // 2 - len(text) // 2, text, curses.A_BOLD)
        self.win_blinker.refresh()

    def _sleep_with_exit(self, duration: float) -> bool:
        """
        Función de ayuda (KISS) que duerme por 'duration'
        mientras comprueba si el usuario presiona 'q' para salir.
        """
        steps = int(duration / 0.05) # Revisar 20 veces por segundo
        if steps == 0: steps = 1

        for _ in range(steps):
            if self.stdscr.getch() == ord('q'):
                return True # Salida solicitada
            time.sleep(duration / steps)
        return False # No se solicitó salida

    def play_text(self, text: str, translator: MorseTranslator):
        """
        Punto de entrada principal para reproducir un texto.

        Args:
            text (str): El texto a reproducir.
            translator (MorseTranslator): El objeto traductor (DIP).
        """

        # 1. Mostrar el texto estático en la ventana superior
        static_morse = translator.get_static_display(text)
        # Manejar texto largo
        max_w = self.win_static.getmaxyx()[1] - 4
        self.win_static.addstr(3, 4, static_morse[:max_w])
        if len(static_morse) > max_w:
            self.win_static.addstr(4, 4, static_morse[max_w:max_w*2])
        self.win_static.refresh()

        # 2. Iterar y reproducir (lógica KISS)
        for i, char in enumerate(text.upper()):
            morse_char = translator.MORSE_CODE_DICT.get(char)

            if morse_char == '/': # Espacio de palabra
                # Un espacio de palabra es 7T. Ya esperamos 3T (letra)
                # así que solo esperamos 4T adicionales.
                self._update_blinker("  /  ")
                if self._sleep_with_exit(self.T_WORD_SPACE - self.T_LETTER_SPACE): return

            elif morse_char: # Es una letra
                for j, symbol in enumerate(morse_char):
                    # A. Mostrar el SÍMBOLO (*)
                    self._update_blinker("*")
                    duration = self.T_DOT if symbol == '.' else self.T_DASH
                    if self._sleep_with_exit(duration): return

                    # B. Mostrar el ESPACIO entre símbolos (1T)
                    self._update_blinker(" ")
                    if self._sleep_with_exit(self.T_SYMBOL_SPACE): return

                # C. Mostrar el ESPACIO entre letras (3T)
                # (Ya esperamos 1T, así que esperamos 2T más)
                if i < len(text) - 1 and text[i+1] != ' ':
                    if self._sleep_with_exit(self.T_LETTER_SPACE - self.T_SYMBOL_SPACE): return

        # 3. Finalizado
        self._update_blinker("[Transmisión completada. 'q' para salir.]")
        while self.stdscr.getch() != ord('q'):
            time.sleep(0.1)


# --- Módulo 3: Aplicación Principal (OCP/DIP) ---
# Esta función une los componentes.

def main_app(stdscr, text_to_play: str):
    """
    Función principal envuelta por curses.wrapper.
    Configura e inyecta las dependencias.
    """

    # --- Configuración (OCP) ---
    # ¡Ajusta estos valores!
    DOT_DURATION = 0.25  # Duración del punto (solicitado: 0.25s)
    DASH_DURATION = 1.0  # Duración de la raya (solicitado: 1.0s)

    # --- Inyección de Dependencias (DIP) ---
    translator = MorseTranslator()
    player = MorsePlayer(stdscr, DOT_DURATION, DASH_DURATION)

    try:
        player.play_text(text_to_play, translator)
    except curses.error:
        # Maneja errores (p.ej., si se redimensiona la terminal)
        pass
    except KeyboardInterrupt:
        pass


# --- Punto de Entrada ---
if __name__ == "__main__":

    # (KISS) Mantenemos el ejemplo simple y contenido
    texto_ejemplo = "SOS HABI TEST"

    if len(sys.argv) > 1:
        texto_ejemplo = " ".join(sys.argv[1:])

    print(f"Iniciando transmisor Morse para: '{texto_ejemplo}'")
    print("La aplicación tomará control de la terminal...")
    print("Presiona 'q' en cualquier momento para salir.")
    time.sleep(2.5)

    # curses.wrapper maneja la inicialización y limpieza
    # de la terminal de forma segura.
    curses.wrapper(main_app, texto_ejemplo)

    print("Transmisor detenido. Terminal restaurada.")
