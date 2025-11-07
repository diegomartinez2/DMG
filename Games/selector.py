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
"""Devuelve lista de posiciones. [] = abortar, no vacío = confirmar cadena."""


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
            elif key == ord(' '):  # <--- NUEVO: SPACE confirma la cadena
                if len(self._selected) >= 2:
                    result = self._selected.copy()
                    self._selected.clear()
                    self._draw()
                    return result
                # Si menos de 2, ignora (o podrías mostrar mensaje)
            elif key == ord('q'):
                return []                # abortar
            self._draw()
