#!/usr/bin/env python3
# laundry-demonology-simulator.py
# Simulador oficial de Computational Demonology - Clasificación: LEVEL 4 / NECRONOMICON CLEARANCE REQUIRED
# Autor: Agent H. (Habi, La Hada del Bikini Azul) - No ejecutar en viernes

import curses
import time
import random
import string
import threading
from itertools import cycle

OPTIONS = {
    "1": ("Protection", "Rosencraft toroidal sequence", "Slow ectastic entity advance"),
    "2": ("Summon", "Level-7 imaginary prime series", "Attract mesonic-energy entities"),
    "3": ("Cast", "nth-cardinality polynomial set", "Disturb JKR boson fields"),
    "4": ("Ftaghu", "Quintic reciprocal convergence", "Tune Planck constant (DANGEROUS)"),
    "5": ("Expel", "Golberg entirion matrices", "Repel boundary anomalies"),
    "Q": ("Exit", "Terminate session", "Please"),
}

def fake_progress(stdscr, y, x, length=40):
    chars = cycle("░▒▓█▓▒░")
    for _ in range(random.randint(20, 80)):
        bar = "".join(next(chars) for _ in range(random.randint(1, length)))
        stdscr.addstr(y, x, f"[{bar.ljust(length)}] {random.randint(1, 100)}%")
        stdscr.refresh()
        time.sleep(0.07)

def scrolling_garbage(stdscr, row):
    while threading.main_thread().is_alive():
        line = ''.join(random.choice("01∞ΣλΔψΩΘΞΦΨϑϒϖℏ∂∇∮∫∰∯∮∬∭∮∯") for _ in range(100))
        stdscr.addstr(row, 0, line[:curses.COLS-1].ljust(curses.COLS-1), curses.A_DIM)
        stdscr.refresh()
        time.sleep(0.15)
        row += 1
        if row >= curses.LINES - 10:
            row = 10

def main(stdscr):
    curses.curs_set(0)
    curses.start_color()
    curses.init_pair(1, curses.COLOR_CYAN, curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_GREEN, curses.COLOR_BLACK)
    curses.init_pair(3, curses.COLOR_RED, curses.COLOR_BLACK)
    curses.init_pair(4, curses.COLOR_YELLOW, curses.COLOR_BLACK)
    stdscr.timeout(100)

    # Fondo de símbolos flotantes
    t = threading.Thread(target=scrolling_garbage, args=(stdscr, 10), daemon=True)
    t.start()

    while True:
        stdscr.clear()
        h, w = stdscr.getmaxyx()

        title = "═══ LAUNDRY COMPUTATIONAL DEMONOLOGY TERMINAL v9.3 (CLASS IV) ═══"
        stdscr.addstr(1, max(0, (w-len(title))//2), title, curses.A_BOLD | curses.color_pair(1))

        stdscr.addstr(3, 2, "CASE NIGHTMARE GREEN status: ", curses.color_pair(4))
        stdscr.addstr("YELLOW/AMBRE (nominal)" if random.random() > 0.2 else "RED/VERMILLION (oh no)", curses.A_BOLD)

        stdscr.addstr(5, 2, "Select necromantic operation:", curses.A_BOLD)
        for key, (name, algo, purpose) in OPTIONS.items():
            if key != "Q":
                color = curses.color_pair(2) if key in "1235" else curses.color_pair(3)
                stdscr.addstr(6 + int(key)-1, 4, f"[{key}] {name:<12} → {algo}", color)
                stdscr.addstr(7 + int(key)-1, 8, f"    └─ {purpose}", curses.A_DIM)

        stdscr.addstr(12, 4, "[Q] Exit / Burn laptop / Run", curses.color_pair(1))

        stdscr.addstr(h-3, 2, "WARNING: Unauthorized use may attract attention of Blue Hades, K-Selected entities or HR.", curses.color_pair(3))
        stdscr.addstr(h-2, 2, f"User: {random.choice(['boble', 'mo', 'andrews', 'angleton', 'guest'])} | Clearance: {random.choice(['LEVEL 2', 'LEVEL 3', 'LEVEL 4', 'EXTERNAL'])}", curses.A_DIM)

        stdscr.refresh()

        key = stdscr.getch()
        if chr(key).upper() == 'Q':
            break
        elif chr(key) in OPTIONS and chr(key) != 'Q':
            operation = OPTIONS[chr(key).upper()]
            stdscr.clear()
            stdscr.addstr(2, 2, f"Iniciando {operation[0]}...", curses.A_BOLD | curses.color_pair(3))
            stdscr.addstr(4, 2, f"Ejecutando: {operation[1]}")
            stdscr.addstr(5, 2, f"Propósito: {operation[2]}")
            stdscr.refresh()
            time.sleep(1.5)

            # Simulación de cálculo pesado
            for i in range(8, 18):
                stdscr.addstr(i, 4, f"Phase {i-7}: Resolving {random.choice(['hypergeometric', 'transfinite', 'non-euclidean', 'necroflux'])} tensors...")
                fake_progress(stdscr, i, 50)
                time.sleep(0.3)

            # Resultado dramático
            result = random.choice([
                "Sequence converged. Local reality stable (±0.0003%).",
                "Minor thaumic backwash detected. Your nose is bleeding.",
                "Entity contact established. It says 'ph'nglui mglw'nafh'. Rebooting.",
                "Planck constant detuned 0.0000000004%. Schrödinger's cat just unionized.",
                "JKR field inversion successful. Tomorrow's lottery numbers: 4 8 15 16 23 42",
                "ERROR: Class-4 breach contained. Angleton wants a word.",
                "Success. The stars are slightly less right today."
            ])
            stdscr.addstr(20, 2, result, curses.A_BOLD | curses.color_pair(2) if "Success" in result else curses.color_pair(3))
            stdscr.addstr(22, 2, "Press any key to continue...")
            stdscr.nodelay(False)
            stdscr.getch()
            stdscr.nodelay(True)

curses.wrapper(main)
