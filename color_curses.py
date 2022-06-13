#colores
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
