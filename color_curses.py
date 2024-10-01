#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
#
#  Copyright 2024 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
##colores
# ---------------------------
# Importación de los módulos
# ---------------------------
import curses

def main(stdscr):
    stdscr.clear()
    if curses.has_colors():
        for i in range(1, curses.COLORS,1):
            curses.init_pair(i, i, curses.COLOR_BLACK)

            stdscr.addstr("COLOR %d! " % i, curses.color_pair(i))
            stdscr.addstr("BOLD! ",         curses.color_pair(i) | curses.A_BOLD)
            stdscr.addstr("STANDOUT! ",     curses.color_pair(i) | curses.A_STANDOUT)
            stdscr.addstr("UNDERLINE! ",    curses.color_pair(i) | curses.A_UNDERLINE)
            stdscr.addstr("BLINK! ",        curses.color_pair(i) | curses.A_BLINK)
            stdscr.addstr("DIM! ",          curses.color_pair(i) | curses.A_DIM)
            stdscr.addstr("REVERSE! ",      curses.color_pair(i) | curses.A_REVERSE)
    stdscr.refresh()
    stdscr.getch()
if __name__ == '__main__':
    print ("init...")
    curses.wrapper(main)
