#!/usr/bin/env python
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
