#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ncurses_clock.py
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
#
#TUI example (making a clock is a way to see how this works)
import curses
import time
import threading


class Clock(threading.Thread):
    """ Clock curses string class. Updates every second. Easy to install """

    def __init__(self, stdscr, show_seconds=True):
        """ Create the clock """
        super(Clock, self).__init__()
        if show_seconds:
            self._target=self.update_seconds
        else:
            self._target=self.blink_colon
        self.daemon = True
        self.stdscr = stdscr
        self.start()

    def update_seconds(self):
        """ If seconds are showing, update the clock each second """
        while 1:
            self.stdscr.addstr(12, 12, time.strftime("%a, %d %b %Y %H:%M:%S"))
            self.stdscr.refresh()
            time.sleep(1)

    def blink_colon(self):
        """ If seconds are not showing, blink the colon each second """
        while 1:
            if int(time.time()) % 2 != 0:
                self.stdscr.addstr(14, 12, time.strftime("%a, %d %b %Y %H:%M"))
            else:
                self.stdscr.addstr(14, 12, time.strftime("%a, %d %b %Y %H %M"))
            stdscr.refresh()
            time.sleep(1)


def run(stdscr):

    stdscr.addstr( 1, 0, "This is sample text\n\n")
    stdscr.addstr(18, 0, "This is more sample text\n\n")
    clock1 = Clock(stdscr)
    clock2 = Clock(stdscr, show_seconds=False)

    # End with any key

    while 1:
        event = stdscr.getch()
        break


if __name__=="__main__":
    stdscr = curses.initscr()
    curses.wrapper(run)
