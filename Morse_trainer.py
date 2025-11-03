import curses
import time

# --- Módulo 1: Traductor (SRP) ---
# Esta clase (reutilizada y mejorada) maneja AMBAS traducciones.
# SRP: Su única responsabilidad es conocer el alfabeto Morse y
# proporcionar métodos para traducir en ambas direcciones.

class MorseTranslator:
#    """
#    Traduce cadenas de texto de/hacia el código Morse Internacional.
#    """

    # Diccionario principal (Texto -> Morse)
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
        ' ': '/' # Espacio entre palabras
    }

    def __init__(self):
        # (KISS) Creamos el diccionario inverso una sola vez.
        self._generate_reverse_dict()

    def _generate_reverse_dict(self):
#        """
#        Función interna (SRP) para crear el diccionario de
#        traducción inversa (Morse -> Texto).
#        """
        self.REVERSE_MORSE_DICT = {}
        for char, code in self.MORSE_CODE_DICT.items():
            self.REVERSE_MORSE_DICT[code] = char

    def get_static_display(self, text: str) -> str:
#        """
#        Traduce texto a una cadena estática de Morse (para el script anterior).
#        """
        morse_text = ""
        for char in text.upper():
            morse_char = self.MORSE_CODE_DICT.get(char, ' ')
            if morse_char == '/':
                morse_text += " /  "
            elif morse_char != ' ':
                morse_text += morse_char + " "
        return morse_text.strip()

    def translate_morse_char(self, morse_code: str) -> str:

#"""
#        Traduce una cadena de morse (ej. '.-') a su letra (ej. 'A').
#
#        Args:
#            morse_code (str): La cadena de puntos y rayas.
#
#        Returns:
#            str: La letra correspondiente, o '?' si no se encuentra.
#"""

        return self.REVERSE_MORSE_DICT.get(morse_code, '?')


# --- Módulo 2: Entrenador Interactivo (SRP) ---
# Esta clase maneja toda la E/S de la terminal (curses) y el estado de la UI.
# No sabe nada del alfabeto, solo usa el traductor.

class MorseTrainer:
#    """
#    Maneja la interfaz de Curses para la práctica interactiva de Morse.
#    Aplica SRP: Su única responsabilidad es la UI y el manejo de estado.
#    """

    def __init__(self, stdscr, translator: MorseTranslator):
#        """
#        Inicializa el entrenador.
#
#        Args:
#            stdscr: La ventana principal de curses.
#            translator (MorseTranslator): El objeto traductor (DIP).
#        """
        self.stdscr = stdscr
        self.translator = translator # Inyección de dependencias

        # --- Estado de la aplicación (KISS) ---
        self.current_morse_buffer = "" # El '.-' que el usuario está escribiendo
        self.current_word = ""         # La palabra que se está formando
        self.full_text = ""            # El texto completo fijado

        # --- Configuración de Curses ---
        curses.curs_set(0)  # Ocultar el cursor
        self.stdscr.nodelay(True) # No bloquear en getch()
        self.stdscr.timeout(50)   # Revisar entrada rápidamente

    def _draw_ui(self):
#        """
#        (SRP) Única función responsable de dibujar la pantalla completa.
#        (KISS) Mantiene la lógica de dibujo en un solo lugar.
#        """
        self.stdscr.clear()
        max_y, max_x = self.stdscr.getmaxyx()

        # 1. Título e Instrucciones
        self.stdscr.addstr(1, 2, "Entrenador de Código Morse Interactivo")
        self.stdscr.addstr(3, 2, "Usa '.' y '-' para formar letras. Pulsa [Espacio] para fijar.")
        self.stdscr.addstr(4, 2, "Pulsa [Espacio] dos veces para un espacio. Pulsa 'q' para salir.")

        # 2. El texto completo que ya se ha fijado
        self.stdscr.addstr(7, 4, "Texto: ", curses.A_BOLD)
        self.stdscr.addstr(self.full_text)

        # 3. La palabra que se está formando actualmente
        self.stdscr.addstr(8, 4, "Palabra: ", curses.A_BOLD)
        self.stdscr.addstr(self.current_word)

        # 4. La entrada dinámica (el núcleo de tu solicitud)
        current_letter = self.translator.translate_morse_char(self.current_morse_buffer)
        if not self.current_morse_buffer:
            current_letter = "_" # Placeholder

        morse_line = f"Búfer Morse: {self.current_morse_buffer}"
        letter_line = f"-> Letra: {current_letter}"

        self.stdscr.addstr(10, 6, morse_line)
        self.stdscr.addstr(10, 6 + len(morse_line) + 2, letter_line, curses.A_REVERSE)

        # 5. Línea de estado
        self.stdscr.addstr(max_y - 2, 2, "Estado: Listo")
        self.stdscr.refresh()

    def _handle_input(self, key: int):
#        """(SRP) Única función responsable de manejar la lógica de entrada."""

        # A. Añadir a búfer (Punto o Raya)
        if key == ord('.') or key == ord('-'):
            self.current_morse_buffer += chr(key)

        # B. Fijar letra (Espacio)
        elif key == ord(' '):
            if self.current_morse_buffer:
                # Había un código en el búfer, fijar la letra
                letter = self.translator.translate_morse_char(self.current_morse_buffer)
                if letter != '?':
                    self.current_word += letter
            else:
                # El búfer estaba vacío (doble espacio)
                # Fijar la palabra actual al texto y añadir un espacio
                self.full_text += self.current_word + " "
                self.current_word = ""

            # Limpiar el búfer para la siguiente letra
            self.current_morse_buffer = ""

    def run_loop(self):
#        """Bucle principal de la aplicación."""
        while True:
            # 1. Dibujar la UI en cada fotograma
            self._draw_ui()

            # 2. Obtener entrada
            try:
                key = self.stdscr.getch()
            except curses.error:
                time.sleep(0.05)
                continue

            # 3. Salir
            if key == ord('q'):
                break

            # 4. Procesar la entrada
            self._handle_input(key)

# --- Módulo 3: Aplicación Principal (OCP/DIP) ---
def main_app(stdscr):
#    """
#    Función principal envuelta por curses.wrapper.
#    Configura e inyecta las dependencias.
#    """
    try:
        # --- Inyección de Dependencias (DIP) ---
        translator = MorseTranslator()
        trainer = MorseTrainer(stdscr, translator)

        # Iniciar el bucle
        trainer.run_loop()

    except KeyboardInterrupt:
        pass
    except curses.error:
        # Maneja errores (p.ej., si se redimensiona la terminal)
        print("Error de Curses. Redimensiona la terminal e inténtalo de nuevo.")


# --- Punto de Entrada ---
if __name__ == "__main__":
    print("Iniciando entrenador Morse...")
    time.sleep(1)
    # curses.wrapper maneja la inicialización y limpieza
    # de la terminal de forma segura.
    curses.wrapper(main_app)
    print("Entrenador detenido. Terminal restaurada.")

#```
#
# ### Cómo usar el script
#
# 1.  **Guardar:** Guarda el código como `morse_trainer.py`.
# 2.  **Ejecutar:**
#    ```bash
#    python3 morse_trainer.py
#'''
