#!/usr/bin/env python3
# laundry-demonology-simulator.py
# Versión corregida – funciona en Linux, macOS y Windows (con windows-curses si hace falta)

import curses
import time
import random
import threading
from itertools import cycle
import sys

OPTIONS = {
    "1": ("Protection", "Rosencraft toroidal sequence", "Slow ectastic entity advance"),
    "2": ("Summon", "Level-7 imaginary prime series", "Attract mesonic-energy entities"),
    "3": ("Cast", "nth-cardinality polynomial set", "Disturb JKR boson fields"),
    "4": ("Ftaghu", "Quintic reciprocal convergence", "Tune Planck constant (DANGEROUS)"),
    "5": ("Expel", "Golberg entirion matrices", "Repel boundary anomalies"),
}

def fake_progress(stdscr, y, x, length=40):
    chars = cycle("░▒▓█▓▒░")
    for _ in range(random.randint(20, 70)):
        bar = "".join(next(chars) for _ in range(random.randint(1, length)))
        stdscr.addstr(y, x, f"[{bar.ljust(length)}] {random.randint(1, 100)}%")
        stdscr.refresh()
        time.sleep(0.07)

def scrolling_garbage(stop_event, stdscr, start_row):
    while not stop_event.is_set():
        line = ''.join(random.choice("01∞ΣλΔψΩΘΞΦΨϑϒϖℏ∂∇∮∫∰∯∮∬∭∮∯∫∲∳") for _ in range(200))
        try:
            stdscr.addstr(start_row, 0, line[:curses.COLS-1].ljust(curses.COLS-1), curses.A_DIM)
        except curses.error:
            pass
        stdscr.refresh()
        time.sleep(0.18)
        start_row += 1
        if start_row >= curses.LINES - 12:
            start_row = 10

def main(stdscr):
    curses.curs_set(0)           # Oculta el cursor
    curses.start_color()
    curses.init_pair(1, curses.COLOR_CYAN, curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_GREEN, curses.COLOR_BLACK)
    curses.init_pair(3, curses.COLOR_RED, curses.COLOR_BLACK)
    curses.init_pair(4, curses.COLOR_YELLOW, curses.COLOR_BLACK)
    stdscr.timeout(100)          # getch() no bloqueante

    stop_event = threading.Event()
    garbage_thread = threading.Thread(target=scrolling_garbage, args=(stop_event, stdscr, 10), daemon=True)
    garbage_thread.start()

    while True:
        stdscr.clear()
        h, w = stdscr.getmaxyx()

        title = " LAUNDRY COMPUTATIONAL DEMONOLOGY TERMINAL v9.3 (CLASS IV) "
        stdscr.addstr(1, max(0, (w-len(title))//2), title, curses.A_BOLD | curses.color_pair(1))

        stdscr.addstr(3, 2, "CASE NIGHTMARE GREEN status: ", curses.color_pair(4))
        status = "YELLOW/AMBER" if random.random() > 0.15 else "RED/VERMILLION – EVACUATE"
        stdscr.addstr(3, 32, status, curses.A_BOLD | (curses.color_pair(3) if "RED" in status else curses.color_pair(2)))

        stdscr.addstr(5, 2, "Select operation:", curses.A_BOLD)
        for i, (key, (name, algo, purpose)) in enumerate(OPTIONS.items(), 7):
            color = curses.color_pair(2) if key != "4" else curses.color_pair(3)
            stdscr.addstr(i, 4, f"[{key}] {name:<12} → {algo}", color)
            stdscr.addstr(i+1, 8, f"    └─ {purpose}", curses.A_DIM)

        stdscr.addstr(h-5, 4, "[Q] Exit and wipe RAM cache", curses.color_pair(1))
        stdscr.addstr(h-3, 2, "WARNING: This terminal may be monitored by Residual Human Resources.", curses.color_pair(3))
        stdscr.addstr(h-2, 2, f"User: {random.choice(['boble','mo','angleton','harry','guest'])} | Clearance: LEVEL {random.randint(2,4)}", curses.A_DIM)

        stdscr.refresh()

        key = stdscr.getch()

        # Protección contra teclas raras
        if key == -1:                    # timeout, nada pulsado
            continue
        if key in (curses.KEY_RESIZE,):  # redimensionar ventana
            continue
        if key >= 0x110000:              # valores fuera de rango
            continue

        try:
            char = chr(key).upper()
        except ValueError:
            continue

        if char == 'Q':
            stop_event.set()
            return

        if char in OPTIONS:
            operation = OPTIONS[char]
            stdscr.clear()
            stdscr.addstr(2, 2, f"INITIALIZING {operation[0].upper()}...", curses.A_BOLD | curses.color_pair(3))
            stdscr.addstr(4, 2, f"Algorithm: {operation[1]}")
            stdscr.addstr(5, 2, f"Purpose:   {operation[2]}")
            stdscr.refresh()
            time.sleep(1.8)

            for phase in range(1, 9):
                stdscr.addstr(7+phase, 4, f"Phase {phase}: Computing {random.choice(['transfinite','hyperelliptic','necromantic','non-commutative'])} coefficients...")
                fake_progress(stdscr, 7+phase, 50)
                time.sleep(0.25)

            results = [
                "Local reality gradient stabilized (±0.0002%). Well done, agent.",
                "Minor thaumic leakage contained. Your soul is slightly singed.",
                "Entity acknowledged presence. It is… curious.",
                "Planck constant jitter: 0.0000000007%. Cats are both alive and annoyed.",
                "JKR field successfully inverted. You just won tomorrow’s lottery (don’t tell anyone).",
                "Class-3 breach averted. Angleton sends his regards.",
                "The stars are almost right… wait, no, false alarm. Phew.",
                "Fhtagn-level resonance detected. Rebooting in 3… 2…",
            ]
            result = random.choice(results)
            color = curses.color_pair(2) if any(x in result for x in ["stabilized","contained","inverted","averted"]) else curses.color_pair(3)
            stdscr.addstr(20, 2, result, curses.A_BOLD | color)
            stdscr.addstr(22, 2, "Press any key to return to menu...", curses.color_pair(1))
            stdscr.nodelay(False)
            stdscr.getch()
            stdscr.nodelay(True)

    stop_event.set()

if __name__ == "__main__":
    try:
        curses.wrapper(main)
    except KeyboardInterrupt:
        print("\nSession terminated by user. Remember to wipe the SSD with a Class-3 exorcism.")
        sys.exit(0)
