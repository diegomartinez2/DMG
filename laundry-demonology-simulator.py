#!/usr/bin/env python3
# laundry-demonology-terminal.py
# Versión FINAL – FUNCIONA EN TODOS LOS SISTEMAS

import curses
import time
import random
import threading
import sys

OPTIONS = {
    "1": ("Protection", "Rosencraft toroidal sequence", "Slow ectastic entity advance"),
    "2": ("Summon", "Level-7 imaginary prime series", "Attract mesonic-energy entities"),
    "3": ("Cast", "nth-cardinality polynomial set", "Disturb JKR boson fields"),
    "4": ("Ftaghu", "Quintic reciprocal convergence", "Wake The Sleeper – DO NOT"),
    "5": ("Expel", "Golberg entirion matrices", "Repel boundary anomalies"),
}

WARNINGS = [
    "WARNING: Unauthorized use may violate the Benthic Treaty",
    "ALERT: Minor thaumic flux in sub-basement 12",
    "NOTICE: Form 27B/6 (Soul Transfer) required in triplicate",
    "CASE NIGHTMARE GREEN readiness level: 71% and dropping",
    "The Auditors are watching. Smile.",
    "Residual Human Resources reminds you: lunch is a privilege",
]

# Variable global para detener el hilo de símbolos
stop_symbols = threading.Event()

def draw_symbols(win):
    symbols = "ΔΘΞΛΨΩ∇∂∮∯∰∭∫∬ℏλ∞∑∝∴∵§¶†‡"
    while not stop_symbols.is_set():
        try:
            h, w = win.getmaxyx()
            for y in range(h):
                line = ''.join(random.choice(symbols) for _ in range(w))
                win.addstr(y, 0, line[:w], curses.A_DIM)
            win.refresh()
        except:
            pass
        time.sleep(0.3)

def print_log(win, lines):
    win.erase()
    h, w = win.getmaxyx()
    for i, line in enumerate(lines[-h+1:]):
        if i >= h - 1:
            break
        color = curses.color_pair(2)
        if any(w in line.upper() for w in ["ERROR", "BREACH", "SLEEPER", "NYARLATHOTEP"]):
            color = curses.color_pair(3)
        elif "SUCCESS" in line.upper():
            color = curses.color_pair(2) | curses.A_BOLD
        win.addstr(i, 0, line[:w-1].ljust(w-1), color)
    win.refresh()

def main(stdscr):
    global stop_symbols
    curses.curs_set(0)
    curses.start_color()
    curses.init_pair(1, curses.COLOR_CYAN,    curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_GREEN,   curses.COLOR_BLACK)
    curses.init_pair(3, curses.COLOR_RED,     curses.COLOR_BLACK)
    curses.init_pair(4, curses.COLOR_YELLOW,  curses.COLOR_BLACK)

    stdscr.nodelay(True)

    while True:
        h, w = stdscr.getmaxyx()

        # Ventanas
        header   = stdscr.derwin(1, w, 0, 0)
        menu     = stdscr.derwin(14, w, 1, 0)
        sep1     = stdscr.derwin(1, w, 15, 0)
        logwin   = stdscr.derwin(h-19, w, 16, 0)
        sep2     = stdscr.derwin(1, w, h-3, 0)
        warn     = stdscr.derwin(1, w, h-2, 0)
        status   = stdscr.derwin(1, w, h-1, 0)

        # Header
        header.addstr(0, 0, "▓" * w, curses.A_REVERSE | curses.color_pair(1))
        header.refresh()

        # Menú
        menu.erase()
        menu.addstr(0, 2, "LAUNDRY – COMPUTATIONAL DEMONOLOGY TERMINAL v10", curses.A_BOLD | curses.color_pair(1))
        case_status = random.choices(["YELLOW/AMBER", "RED/VERMILLION"], weights=[93,7])[0]
        menu.addstr(2, 2, "CASE NIGHTMARE GREEN: ", curses.color_pair(4))
        menu.addstr(2, 24, case_status, curses.A_BOLD | (curses.color_pair(3) if "RED" in case_status else curses.color_pair(2)))

        menu.addstr(4, 2, "Select operation:", curses.A_BOLD)
        for i, (k, (name, algo, _)) in enumerate(OPTIONS.items(), 6):
            color = curses.color_pair(3) | curses.A_BOLD if k == "4" else curses.color_pair(2)
            menu.addstr(i, 4, f"[{k}] {name:<12} → {algo}", color)

        menu.addstr(12, 4, "[Q] Exit and burn evidence", curses.color_pair(1))
        menu.refresh()

        # Separadores
        sep1.addstr(0, 0, "─" * w, curses.A_DIM)
        sep2.addstr(0, 0, "─" * w, curses.A_DIM)
        sep1.refresh(); sep2.refresh()

        # Warning y status
        warn.addstr(0, 0, random.choice(WARNINGS), curses.color_pair(4) | curses.A_BOLD)
        warn.refresh()

        user = random.choice(["boble", "mo", "angleton", "harriet", "persephone"])
        status.addstr(0, 0, f"User: {user.upper()} │ Clearance: LEVEL {random.randint(2,4)} │ Entropy: {random.randint(40,98)}%")
        status.refresh()

        # Iniciar hilo de símbolos solo la primera vez
        if not stop_symbols.is_set() and threading.active_count() < 3:
            threading.Thread(target=draw_symbols, args=(logwin,), daemon=True).start()

        print_log(logwin, ["Ready.", "Awaiting command..."])

        # Entrada
        key = stdscr.getch()
        if key == -1:
            continue
        if key in (curses.KEY_RESIZE,):
            continue

        try:
            ch = chr(key).upper()
        except:
            continue

        if ch == 'Q':
            stop_symbols.set()
            break

        if ch in OPTIONS:
            op = OPTIONS[ch]
            log = [
                f"> INITIATING {op[0].upper()}...",
                f"  Algorithm: {op[1]}",
                f"  Purpose: {op[2]}",
                ""
            ]
            print_log(logwin, log)

            for phase in range(1, 9):
                log.append(f"  Phase {phase}/8: {random.choice(['Necromantic recursion','Transfinite tensor collapse','Hyperdimensional inversion'])} in progress...")
                print_log(logwin, log)
                time.sleep(0.35 + random.random() * 0.4)

            result = random.choice([
                "SUCCESS: Local manifold stabilized.",
                "MINOR BREACH: Your soul now belongs 6.66% to an outer god.",
                "ERROR: The Sleeper has opened one eye. Run.",
                "JKR field inverted. You just rewrote Tuesday.",
                "Entity repelled. Angleton owes you a beer.",
                "Planck constant successfully detuned. All physics now optional.",
            ])
            log.append("")
            log.append(result)
            print_log(logwin, log)
            time.sleep(3.5)

    # Salida
    stdscr.clear()
    msg = "TERMINAL SELF-DESTRUCT SEQUENCE INITIATED"
    stdscr.addstr(h//2, (w-len(msg))//2, msg, curses.A_BOLD | curses.color_pair(3))
    stdscr.refresh()
    time.sleep(3)

if __name__ == "__main__":
    try:
        curses.wrapper(main)
    except KeyboardInterrupt:
        stop_symbols.set()
        print("\nEmergency shutdown. Please report to Decontamination.")
        sys.exit(0)
