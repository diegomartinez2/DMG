# -------------------------------------------------
# file: game.py
# -------------------------------------------------
from typing import List, Tuple
from board import IBoard, Board
from selector import ISelector, CursorSelector
from scorer import IScorer, LengthScorer
import curses


class Game:
    """
    Orquestador del juego.
    - Depende únicamente de abstracciones (IBoard, ISelector, IScorer).
    - Fácil de probar o sustituir cualquier componente.
    """

    def __init__(self, board: IBoard, selector: ISelector, scorer: IScorer):
        self._board = board
        self._selector = selector
        self._scorer = scorer
        self._total = 0

    # ------------------------------------------------------------------
    def _print_result(self, chain: List[Tuple[int, int]], points: int):
        rows, cols = self._board.size()
        stdscr = self._selector._stdscr if hasattr(self._selector, "_stdscr") else None
        if stdscr:
            y = rows + 3
            stdscr.addstr(y, 0, f"Cadena: {' → '.join(self._board.get(r, c) for r, c in chain)}")
            stdscr.addstr(y + 1, 0, f"Puntos: {points}   Total acumulado: {self._total}")
            stdscr.addstr(y + 3, 0, "Presiona cualquier tecla para continuar...")
            stdscr.refresh()
            stdscr.getch()

    # ------------------------------------------------------------------
    def run(self):
        stdscr = curses.initscr()
        curses.noecho()
        curses.cbreak()
        stdscr.keypad(True)

        # Inyectamos stdscr al selector concreto
        selector = CursorSelector(self._board, stdscr)
        self._selector = selector

        try:
            while True:
                chain = selector.run()
                if not chain:                 # abortado con 'q'
                    break
                points = self._scorer.score(chain)
                self._total += points
                self._print_result(chain, points)
        finally:
            curses.endwin()

        print(f"\n¡Juego terminado! Puntuación final: {self._total}")
