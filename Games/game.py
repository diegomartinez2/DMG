# -------------------------------------------------
# file: game.py  (versión actualizada)
# -------------------------------------------------
from typing import List, Tuple
from board import IBoard
from selector import ISelector, CursorSelector
from scorer import IScorer, LengthScorer
import curses


class Game:
    def __init__(self, board: IBoard, selector: ISelector, scorer: IScorer):
        self._board = board
        self._selector = selector
        self._scorer = scorer
        self._total = 0

    def _print_result(self, chain: List[Tuple[int, int]], points: int):
        rows, _ = self._board.size()
        stdscr = self._selector._stdscr if hasattr(self._selector, "_stdscr") else None
        if stdscr:
            y = rows + 3
            letters = ' → '.join(self._board.get(r, c) for r, c in chain)
            stdscr.addstr(y, 0, f"Cadena: {letters}")
            stdscr.addstr(y + 1, 0, f"Puntos: {points}   Total: {self._total}")
            stdscr.addstr(y + 3, 0, "SPACE: nueva cadena  q: salir")
            stdscr.refresh()
            # Espera confirmación antes de limpiar
            while True:
                k = stdscr.getch()
                if k == ord(' '):
                    break
                if k == ord('q'):
                    raise KeyboardInterrupt  # forzar salida

    def run(self):
        stdscr = curses.initscr()
        curses.noecho()
        curses.cbreak()
        stdscr.keypad(True)

        # Inyectamos stdscr
        selector = CursorSelector(self._board, stdscr)
        self._selector = selector

        try:
            while True:
                chain = selector.run()
                if not chain:  # 'q' presionado
                    break
                points = self._scorer.score(chain)
                self._total += points
                self._print_result(chain, points)
                # El tablero ya se limpia en CursorSelector al confirmar
        except KeyboardInterrupt:
            pass
        finally:
            curses.endwin()

        print(f"\nJuego terminado. Puntuación final: {self._total}")
