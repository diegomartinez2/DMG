#!/usr/bin/env python3
# laundry-demonology-terminal.py
# Versión ULTRA-ESTABLE – funciona en cualquier terminal del planeta

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
    "WARNING: Use may violate Benthic Treaty Article IV",
    "ALERT: Thaumic turbulence detected in sub-basement",
    "NOTICE: Form 27B/6 required for soul transfer",
    "CASE NIGHTMARE GREEN readiness: dropping fast",
    "The Auditors are watching. Do not blink.",
    "R.H.R. reminds you: lunch is not guaranteed",
]

stop_symbols = threading.Event()

def draw_symbols(win):
    chars = ".-~:+*#=%&@"
    while not stop_symbols.is_set():
        try:
            h, w = win.getmaxyx()
            if h < 3 or w < 10:
                time.sleep(0.5)
                continue
            for y in range(h):
                line = ''.join(random.choice(chars) for _ in range(w))
                win.addstr(y, 0, line[:w])
            win.refresh()
        except:
            pass
        time.sleep(0.35)

def print_log(win, lines):
    try:
        win.erase()
        h, w = win.getmaxyx()
        start = max(0, len(lines) - (h-1))
        for i, line in enumerate(lines[start:]):
            if i >= h-1:
                break
            text = line[:w-1].ljust(w-1)
            if any(x in line.upper() for x in ["ERROR","BREACH","SLEEPER","NYARLATHOTEP","DO NOT"]):
                color = curses.color_pair(3) | curses.A_BOLD
            elif "SUCCESS" in line.upper():
                color = curses.color_pair(2) | curses.A_BOLD
            else:
                color = curses.A_NORMAL
            win.addstr(i, 0, text, color)
        win.refresh()
    except:
        pass

def main(stdscr):
    global stop_symbols
    curses.curs_set(0)
    stdscr.nodelay(True)

    # Colores seguros
    if curses.has_colors():
        curses.start_color()
        curses.init_pair(1, curses.COLOR_CYAN,    curses.COLOR_BLACK)
        curses.init_pair(2, curses.COLOR_GREEN,   curses.COLOR_BLACK)
        curses.init_pair(3, curses.COLOR_RED,     curses.COLOR_BLACK)
        curses.init_pair(4, curses.COLOR_YELLOW,  curses.COLOR_BLACK)

    log_lines = ["*** LAUNDRY TERMINAL READY ***", "Wards active. Awaiting orders."]

    # Iniciar hilo de símbolos
    threading.Thread(target=draw_symbols, args=(stdscr.derwin(1,1,0,0),), daemon=True).start()

    while True:
        try:
            h, w = stdscr.getmaxyx()
            if h < 20 or w < 60:           # terminal muy pequeño
                stdscr.erase()
                stdscr.addstr(0, 0, "TERMINAL TOO SMALL - RESIZE WINDOW", curses.A_BOLD)
                stdscr.refresh()
                time.sleep(1)
                continue

            # --- Ventanas ---
            header  = stdscr.derwin(1, w, 0, 0)
            menu    = stdscr.derwin(14, w, 1, 0)
            sep1    = stdscr.derwin(1, w, 15, 0)
            logwin  = stdscr.derwin(h-19, w, 16, 0)
            sep2    = stdscr.derwin(1, w, h-3, 0)
            warn    = stdscr.derwin(1, w, h-2, 0)
            status  = stdscr.derwin(1, w, h-1, 0)

            # Header
            header.addstr(0, 0, "=" * (w-1), curses.A_REVERSE | curses.color_pair(1))
            header.refresh()

            # Menú
            menu.erase()
            menu.addstr(0, 2, "LAUNDRY - COMPUTATIONAL DEMONOLOGY v10", curses.A_BOLD | curses.color_pair(1))
            case = random.choices(["YELLOW/AMBER", "RED/VERMILLION"], weights=[94,6])[0]
            menu.addstr(2, 2, "CASE NIGHTMARE GREEN: " + case, curses.color_pair(4) if "YELLOW" in case else curses.color_pair(3)|curses.A_BOLD)
            menu.addstr(4, 2, "Select operation:", curses.A_BOLD)
            for i, (k, (name, algo, _)) in enumerate(OPTIONS.items(), 6):
                col = curses.color_pair(3)|curses.A_BOLD if k == "4" else curses.color_pair(2)
                menu.addstr(i, 4, f"[{k}] {name:<12} -> {algo}", col)
            menu.addstr(12, 4, "[Q] Exit and destroy evidence", curses.color_pair(1))
            menu.refresh()

            # Separadores
            sep1.addstr(0, 0, "-" * (w-1), curses.A_DIM)
            sep2.addstr(0, 0, "-" * (w-1), curses.A_DIM)
            sep1.refresh(); sep2.refresh()

            # Warning y status
            warn.addstr(0, 0, random.choice(WARNINGS), curses.color_pair(4)|curses.A_BOLD)
            warn.refresh()
            user = random.choice(["boble","mo","angleton","harriet","persephone"])
            status.addstr(0, 0, f"User: {user.upper()} | Clearance: LEVEL {random.randint(2,4)} | Entropy: {random.randint(40,99)}%")
            status.refresh()

            # Símbolos de fondo solo en la zona de log
            threading.Thread(target=draw_symbols, args=(logwin,), daemon=True).start()

            print_log(logwin, log_lines)

            # Entrada
            key = stdscr.getch()
            if key == -1:
                continue

            try:
                ch = chr(key).upper()
            except:
                continue

            if ch == "Q":
                stop_symbols.set()
                break

            if ch in OPTIONS:
                op = OPTIONS[ch]
                log_lines.append("")
                log_lines.append(f">>> EXECUTING {op[0].upper()}")
                log_lines.append(f"    Algorithm: {op[1]}")
                print_log(logwin, log_lines)

                for phase in range(1, 9):
                    log_lines.append(f"    Phase {phase}/8: {random.choice(['Necromantic backprop','Transfinite collapse','Hyperdimensional fold'])}...")
                    print_log(logwin, log_lines)
                    time.sleep(0.4 + random.random()*0.4)

                result = random.choice([
                    "SUCCESS: Reality deviation contained",
                    "MINOR BREACH: Soul ownership contested",
                    "ERROR: The Sleeper stirs. Evacuate.",
                    "JKR inversion complete - Tuesday cancelled",
                    "Entity repelled. Pint owed by Angleton",
                    "Planck constant detuned. Physics optional",
                ])
                log_lines.append("")
                log_lines.append(result)
                print_log(logwin, log_lines)
                time.sleep(3)

        except curses.error:
            pass
        except:
            break

    # Salida final
    stdscr.clear()
    msg = "TERMINAL OFFLINE - EVIDENCE DESTROYED"
    stdscr.addstr(h//2, max(0, (w-len(msg))//2), msg, curses.A_BOLD | curses.color_pair(3))
    stdscr.refresh()
    time.sleep(3)

if __name__ == "__main__":
    try:
        curses.wrapper(main)
    except KeyboardInterrupt:
        stop_symbols.set()
        print("\nEmergency shutdown complete. Report to Debriefing Room 12.")
