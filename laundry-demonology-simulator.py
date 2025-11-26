#!/usr/bin/env python3
# laundry-demonology-terminal-v10.py
# Terminal oficial de Computational Demonology - The Laundry
# Ahora con layout profesional de 4 zonas

import curses
import time
import random
import threading
from itertools import cycle

OPTIONS = {
    "1": ("Protection", "Rosencraft toroidal sequence", "Slow ectastic entity advance"),
    "2": ("Summon", "Level-7 imaginary prime series", "Attract mesonic-energy entities"),
    "3": ("Cast", "nth-cardinality polynomial set", "Disturb JKR boson fields"),
    "4": ("Ftaghu", "Quintic reciprocal convergence", "Observe The Sleeper (EXTREME RISK)"),
    "5": ("Expel", "Golberg entirion matrices", "Repel boundary anomalies"),
}

# Colores
PAIR_TITLE   = 1  # cyan
PAIR_GOOD    = 2  # green
PAIR_DANGER  = 3  # red
PAIR_WARN    = 4  # yellow
PAIR_DIM     = 5  # gris

warnings = [
    "WARNING: Unauthorized observation may violate the Kibble Protocol",
    "ALERT: Minor thaumic turbulence detected in sector 7G",
    "NOTICE: Mahogany Row requests immediate status report",
    "Residual Human Resources reminds you: no soul transfers without Form 27B/6",
    "CASE NIGHTMARE GREEN readiness: 73% (and falling)",
    "The Auditors are performing a surprise inspection in 00:03:14",
]

def scrolling_symbols(win):
    symbols = "ΔΘΞΛΨΩ∇∂∮∯∰∭∫∬ℏλ∞∑∝∴∵"
    while not win.stop_event.is_set():
        line = ''.join(random.choice(symbols) for _ in range(200))
        for row in range(win.height):
            win.addstr(row, 0, line[row:row+win.width-1].ljust(win.width-1), curses.A_DIM)
        win.refresh()
        time.sleep(0.25)

def computation_output(win, text_lines):
    win.clear()
    max_lines = win.height - 1
    for i, line in enumerate(text_lines[-max_lines:]):
        color = curses.color_pair(PAIR_GOOD) if any(x in line.lower() for x in ["success","converged","stable"]) else curses.color_pair(PAIR_DIM)
        if "ERROR" in line or "BREACH" in line:
            color = curses.color_pair(PAIR_DANGER)
        win.addstr(i, 0, line[:win.width-1].ljust(win.width-1), color)
    win.refresh()

def main(stdscr):
    curses.curs_set(0)
    curses.start_color()
    curses.init_pair(PAIR_TITLE,   curses.COLOR_CYAN,    curses.COLOR_BLACK)
    curses.init_pair(PAIR_GOOD,    curses.COLOR_GREEN,   curses.COLOR_BLACK)
    curses.init_pair(PAIR_DANGER, curses.COLOR_RED,     curses.COLOR_BLACK)
    curses.init_pair(PAIR_WARN,    curses.COLOR_YELLOW, curses.COLOR_BLACK)
    curses.init_pair(PAIR_DIM,     curses.COLOR_WHITE,  curses.COLOR_BLACK)
    stdscr.timeout(100)

    h, w = stdscr.getmaxyx()

    # Definición de ventanas
    header_win     = stdscr.derwin(1, w, 0, 0)           # línea 0
    menu_win       = stdscr.derwin(14, w, 1, 0)          # zona mensajes + menú
    separator1     = stdscr.derwin(1, w, 15, 0)
    compute_win    = stdscr.derwin(h-19, w, 16, 0)       # zona cálculo
    separator2     = stdscr.derwin(1, w, h-3, 0)
    warning_win    = stdscr.derwin(1, w, h-2, 0)         # línea de warnings
    status_win     = stdscr.derwin(1, w, h-1, 0)         # línea de estado

    compute_win.stop_event = threading.Event()
    symbol_thread = threading.Thread(target=scrolling_symbols, args=(compute_win,), daemon=True)
    symbol_thread.start()

    log_lines = ["=== SYSTEM ONLINE ===", "All wards nominal.", "Awaiting operator input..."]

    while True:
        # Header con símbolos flotantes
        header_win.clear()
        header_win.addstr(0, 0, "▓" * w, curses.A_REVERSE | curses.color_pair(PAIR_TITLE))
        header_win.refresh()

        # Menú principal
        menu_win.clear()
        menu_win.addstr(0, 2, "LAUNDRY - Computational Demonology Terminal v10.0", curses.A_BOLD | curses.color_pair(PAIR_TITLE))
        menu_win.addstr(2, 2, "CASE NIGHTMARE GREEN status:", curses.color_pair(PAIR_WARN))
        status = random.choices(["YELLOW/AMBER", "RED/VERMILLION – EVACUATE"], weights=[90,10])[0]
        menu_win.addstr(2, 31, status, curses.A_BOLD | (curses.color_pair(PAIR_DANGER) if "RED" in status else curses.color_pair(PAIR_GOOD)))

        menu_win.addstr(4, 2, "Select operation:", curses.A_BOLD)
        for i, (k, (name, algo, _)) in enumerate(OPTIONS.items(), start=6):
            color = curses.color_pair(PAIR_DANGER) if k == "4" else curses.color_pair(PAIR_GOOD)
            menu_win.addstr(i, 4, f"[{k}] {name:<12} → {algo}", color)

        menu_win.addstr(12, 4, "[Q] Exit and initiate thermite self-destruct sequence", curses.color_pair(PAIR_TITLE))
        menu_win.refresh()

        # Separadores
        separator1.addstr(0, 0, "─" * (w-1), curses.color_pair(PAIR_DIM))
        separator2.addstr(0, 0, "─" * (w-1), curses.color_pair(PAIR_DIM))
        separator1.refresh(); separator2.refresh()

        # Warning aleatorio
        warning_win.clear()
        warning_win.addstr(0, 0, random.choice(warnings), curses.color_pair(PAIR_WARN) | curses.A_BOLD)
        warning_win.refresh()

        # Línea de estado
        status_win.clear()
        user = random.choice(["boble", "mo", "angleton", "harriet", "derek"])
        status_win.addstr(0, 0, f"User: {user.upper()} | Clearance: LEVEL {random.randint(2,4)} | Entropy: {random.randint(38,97)}% ", curses.color_pair(PAIR_DIM))
        status_win.refresh()

        # Entrada
        key = stdscr.getch()
        if key in (-1, curses.KEY_RESIZE):
            continue
        try:
            ch = chr(key).upper()
        except:
            continue

        if ch == 'Q':
            compute_win.stop_event.set()
            break

        if ch in OPTIONS:
            op = OPTIONS[ch]
            log_lines.append(f"> Initiating {op[0]}...")
            log_lines.append(f"  Algorithm: {op[1]}")
            computation_output(compute_win, log_lines)

            for phase in range(1, 9):
                log_lines.append(f"  Phase {phase}/8: Resolving {random.choice(['transfinite','non-Euclidean','necromantic','hypergeometric'])} tensors...")
                computation_output(compute_win, log_lines)
                time.sleep(0.35 + random.uniform(0, 0.4))

            result = random.choice([
                "SUCCESS: Local reality gradient stabilized.",
                "MINOR BREACH: Your soul is now 3.7% property of Nyarlathotep.",
                "ERROR: The Sleeper stirs. Run.",
                "JKR field inverted. You just changed tomorrow.",
                "Class-4 entity repelled. Angleton owes you a pint.",
                "Planck constant tuned. Schrödinger's cat just filed a complaint.",
            ])
            log_lines.append(result)
            computation_output(compute_win, log_lines)
            time.sleep(2.5)

    # Mensaje final
    stdscr.clear()
    stdscr.addstr(h//
