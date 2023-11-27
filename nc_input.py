import curses

def calculate(a, b, c):
   return a + b + c

def main(stdscr):
   # Initialize the curses object
   stdscr = curses.initscr()

   # Set the mode of the terminal to non-blocking or cbreak mode
   curses.cbreak()

   # Do not echo keys back to the client
   curses.noecho()

   # Get parameters a, b, and c from the user
   stdscr.addstr(0, 0, "Enter parameter a: ")
   a = int(stdscr.getstr().decode('utf-8'))
   stdscr.addstr(1, 0, "Enter parameter b: ")
   b = int(stdscr.getstr().decode('utf-8'))
   stdscr.addstr(2, 0, "Enter parameter c: ")
   c = int(stdscr.getstr().decode('utf-8'))

   # Call the function with the parameters
   result = calculate(a, b, c)

   # Display the result
   stdscr.addstr(4, 0, "Result: " + str(result))

   # Wait for the user to press a key
   stdscr.getch()

   # Deinitialize the curses mode
   curses.endwin()
if __name__ == "__main__":
      curses.wrapper(main)
