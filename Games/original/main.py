#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  main.py
#
#  Copyright 2026 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
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
