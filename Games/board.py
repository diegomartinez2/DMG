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
