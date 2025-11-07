# -------------------------------------------------
# file: board.py
# -------------------------------------------------
from __future__ import annotations
from typing import Protocol, List, Tuple
import random
import string


class IBoard(Protocol):
    """Interfaz mínima de una cuadrícula de juego."""

    def size(self) -> Tuple[int, int]:
        """(filas, columnas)"""

    def get(self, row: int, col: int) -> str:
        """Letra en la posición dada."""

    def __str__(self) -> str:
        """Representación textual completa."""


class Board:
    """
    Implementación concreta de IBoard.
    - Genera una cuadrícula de letras aleatorias.
    - Permite acceso de solo lectura.
    """

    def __init__(self, rows: int = 6, cols: int = 8, alphabet: str = string.ascii_uppercase):
        if rows < 1 or cols < 1:
            raise ValueError("Dimensiones deben ser positivas")
        self._rows = rows
        self._cols = cols
        self._alphabet = alphabet
        self._grid: List[List[str]] = self._generate_grid()

    # ------------------------------------------------------------------
    # Métodos privados
    # ------------------------------------------------------------------
    def _generate_grid(self) -> List[List[str]]:
        return [
            [random.choice(self._alphabet) for _ in range(self._cols)]
            for _ in range(self._rows)
        ]

    # ------------------------------------------------------------------
    # API pública (IBoard)
    # ------------------------------------------------------------------
    def size(self) -> Tuple[int, int]:
        return self._rows, self._cols

    def get(self, row: int, col: int) -> str:
        if not (0 <= row < self._rows and 0 <= col < self._cols):
            raise IndexError("Posición fuera de rango")
        return self._grid[row][col]

    def __str__(self) -> str:
        lines = []
        for r in range(self._rows):
            lines.append(" ".join(self._grid[r]))
        return "\n".join(lines)


# -------------------------------------------------
# file: selector.py
# -------------------------------------------------
from typing import Protocol, List, Tuple
from board import IBoard
import curses   # manejo de teclado en terminal


class ISelector(Protocol):
    """Interfaz para la selección de celdas contiguas."""

    def run(self) -> List[Tuple[int, int]]:
        """Devuelve lista de (row, col) seleccionados o [] si se cancela."""


class CursorSelector:
    """
    Implementación con curses.
    - Flechas mueven el cursor.
    - ENTER agrega la celda actual a la selección.
    - BACKSPACE elimina la última celda.
    - 'q' aborta.
    """

    def __init__(self, board: IBoard, stdscr):
        self._board = board
        self._stdscr = stdscr
        self._selected: List[Tuple[int, int]] = []
        self._cursor = (0, 0)          # posición actual del cursor

    # ------------------------------------------------------------------
    # Helpers internos
    # ------------------------------------------------------------------
    def _draw(self):
        self._stdscr.clear()
        rows, cols = self._board.size()
        # Dibuja tablero
        for r in range(rows):
            for c in range(cols):
                ch = self._board.get(r, c)
                attr = curses.A_REVERSE if (r, c) in self._selected else 0
                if (r, c) == self._cursor:
                    attr |= curses.A_BOLD
                self._stdscr.addstr(r, c * 2, ch, attr)
        # Instrucciones
        instr = "Flechas: mover  ENTER: marcar  BACKSPACE: deshacer  q: salir"
        self._stdscr.addstr(rows + 1, 0, instr)
        self._stdscr.refresh()

    def _move_cursor(self, dr: int, dc: int):
        r, c = self._cursor
        rows, cols = self._board.size()
        r = (r + dr) % rows
        c = (c + dc) % cols
        self._cursor = (r, c)

    def _toggle_current(self):
        pos = self._cursor
        if pos in self._selected:
            self._selected.remove(pos)
        else:
            # Sólo permite celdas adyacentes a la última seleccionada
            if not self._selected or self._is_adjacent(self._selected[-1], pos):
                self._selected.append(pos)

    def _is_adjacent(self, a: Tuple[int, int], b: Tuple[int, int]) -> bool:
        ra, ca = a
        rb, cb = b
        return max(abs(ra - rb), abs(ca - cb)) == 1

    # ------------------------------------------------------------------
    # API pública (ISelector)
    # ------------------------------------------------------------------
    def run(self) -> List[Tuple[int, int]]:
        curses.curs_set(0)               # ocultar cursor del sistema
        self._draw()

        while True:
            key = self._stdscr.getch()

            if key == curses.KEY_UP:
                self._move_cursor(-1, 0)
            elif key == curses.KEY_DOWN:
                self._move_cursor(1, 0)
            elif key == curses.KEY_LEFT:
                self._move_cursor(0, -1)
            elif key == curses.KEY_RIGHT:
                self._move_cursor(0, 1)
            elif key in (10, 13, curses.KEY_ENTER):   # ENTER
                self._toggle_current()
            elif key == curses.KEY_BACKSPACE or key == 127:
                if self._selected:
                    self._selected.pop()
            elif key == ord('q'):
                return []                # abortar
            self._draw()


# -------------------------------------------------
# file: scorer.py
# -------------------------------------------------
from typing import Protocol, List, Tuple


class IScorer(Protocol):
    def score(self, chain: List[Tuple[int, int]]) -> int:
        """Puntos según longitud de la cadena."""


class LengthScorer:
    """
    Puntuación triangular: n(n-1)/2
    - 2 letras → 1 pt
    - 3 letras → 3 pts
    - 4 letras → 6 pts …
    """

    @staticmethod
    def score(chain: List[Tuple[int, int]]) -> int:
        n = len(chain)
        if n < 2:
            return 0
        return n * (n - 1) // 2


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


# -------------------------------------------------
# file: main.py  (punto de entrada)
# -------------------------------------------------
from game import Game
from board import Board
from selector import CursorSelector
from scorer import LengthScorer


def main():
    # Configuración básica (puede cambiarse sin tocar el resto)
    board = Board(rows=6, cols=10)
    # selector y scorer usan sus implementaciones por defecto
    game = Game(board=board, selector=CursorSelector, scorer=LengthScorer())
    game.run()


if __name__ == "__main__":
    main()
