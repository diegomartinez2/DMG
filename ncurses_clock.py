#!/usr/bin/env python
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
