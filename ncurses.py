import curses

def create_window(rows, cols, y, x, title):
    win = curses.newwin(rows, cols, y, x)

    win.box()
    win.move(0, 1)
    win.addstr(title)
    win.move(1, 1)
    return win

def create_centered_window(rows, cols, title):
    y = (curses.LINES - rows) // 2
    x = (curses.COLS - cols) // 2
    return create_window(rows, cols, y, x, title)

def get_name(win):
    curses.echo()
    result = win.getstr()
    curses.noecho()
    return result

def main(stdscr):
    win = create_centered_window(10, 40, 'Introduzca su nombre')
    name = get_name(win)
    win.erase()

    win = create_centered_window(5, 20, 'Saludo')
    win.addstr('Hola, ')
    win.addstr(name)

    stdscr.refresh()
    win.getch()

if __name__ == '__main__':
    curses.wrapper(main)
