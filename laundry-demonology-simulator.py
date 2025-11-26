#!/usr/bin/env python3
# laundry-demonology-terminal.py
# VERSIÓN FINAL INDESTRUCTIBLE - Funciona en cualquier terminal

import curses
import time
import random
import threading

OPTIONS = {
    "1": ("Protection", "Rosencraft toroidal sequence", "Slow ectastic entity advance"),
    "2": ("Summon", "Level-7 imaginary prime series", "Attract mesonic-energy entities"),
    "3": ("Cast", "nth-cardinality polynomial set", "Disturb JKR boson fields"),
    "4": ("Ftaghu", "Quintic reciprocal convergence", "DO NOT RUN ON FRIDAY"),
    "5": ("Expel", "Golberg entirion matrices", "Repel boundary anomalies"),
}

WARNINGS = [
    "WARNING: Unauthorized use violates Benthic Treaty",
    "ALERT: Thaumic turbulence in sub-basement",
    "Form 27B/6 required for soul transfers",
    "CASE NIGHTMARE GREEN readiness dropping",
    "The Auditors are watching",
    "R.H.R.: coffee is not a food group",
]

stop_symbols = threading.Event()
current_user = "GUEST"
current_clearance = 2

def safe_addstr(win, y, x, text, attr=0):
    try:
        win.addstr(y, x, text[:win.getmaxyx()[1]-x-1], attr)
    except:
        pass

def draw_symbols(win):
    chars = ".-:=+*#%&@"
    while not stop_symbols.is_set():
        try:
            h, w = win.getmaxyx()
            if h < 3 or w < 10:
                time.sleep(0.5); continue
            for y in range(h):
                line = "".join(random.choice(chars) for _ in range(w))
                safe_addstr(win, y, 0, line)
            win.refresh()
        except:
            pass
        time.sleep(0.35)

def print_log(win, lines):
    try:
        win.erase()
        h, w = win.getmaxyx()
        for i, line in enumerate(lines[-(h-1):]):
            if i >= h-1: break
            text = str(line)[:w-1].ljust(w-1)
            if any(x in text.upper() for x in ["ERROR","BREACH","SLEEPER","EVACUATE"]):
                attr = curses.A_BOLD | curses.color_pair(3)
            elif "SUCCESS" in text.upper():
                attr = curses.A_BOLD | curses.color_pair(2)
            else:
                attr = curses.A_NORMAL
            safe_addstr(win, i, 0, text, attr)
        win.refresh()
    except:
        pass

def login_screen(stdscr):
    global current_user, current_clearance
    stdscr.clear()
    curses.curs_set(1)
    stdscr.nodelay(False)

    lines = [
        "",
        "    =========================================",
        "    ==   THE LAUNDRY - OCCULT COMPUTING    ==",
        "    ==   COMPUTATIONAL DEMONOLOGY TERMINAL ==",
        "    =========================================",
        "    CLASSIFIED LEVEL 3 AND ABOVE ONLY",
        "",
        "    Username: ",
        "    Clearance (2-4): ",
    ]
    for i, line in enumerate(lines):
        safe_addstr(stdscr, 3+i, 5, line, curses.A_BOLD | curses.color_pair(1))

    curses.echo()
    user_input = stdscr.getstr(9, 18, 20).decode(errors='ignore').strip()
    cl_input = stdscr.getstr(10, 25, 1).decode(errors='ignore').strip()
    curses.noecho()
    curses.curs_set(0)

    current_user = (user_input or "guest").upper()[:12]
    try:
        current_clearance = int(cl_input)
        if current_clearance not in (2,3,4): raise ValueError
    except:
        current_clearance = 2

    safe_addstr(stdscr, 13, 5, "Authentication successful. Loading...", curses.color_pair(2))
    stdscr.refresh()
    time.sleep(1.5)

def main(stdscr):
    global stop_symbols

    curses.start_color()
    curses.init_pair(1, curses.COLOR_CYAN,    curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_GREEN,   curses.COLOR_BLACK)
    curses.init_pair(3, curses.COLOR_RED,     curses.COLOR_BLACK)
    curses.init_pair(4, curses.COLOR_YELLOW,  curses.COLOR_BLACK)

    stdscr.nodelay(True)
    login_screen(stdscr)

    log_lines = [
        "*** SYSTEM ONLINE ***",
        f"User {current_user} - LEVEL {current_clearance}",
        "Wards active. Ready.",
    ]

    while True:
        try:
            h, w = stdscr.getmaxyx()
            if h < 24 or w < 70:
                stdscr.erase()
                safe_addstr(stdscr, 0, 0, "TERMINAL TOO SMALL - RESIZE (70x24 min)", curses.A_BOLD)
                stdscr.refresh()
                time.sleep(1)
                continue

            # Crear ventanas solo si hay espacio
            header = stdscr.derwin(1, w, 0, 0)
            menu   = stdscr.derwin(14, w, 1, 0)
            sep1   = stdscr.derwin(1, w, 15, 0)
            logwin = stdscr.derwin(h-19, w, 16, 0)
            sep2   = stdscr.derwin(1, w, h-3, 0)
            warn   = stdscr.derwin(1, w, h-2, 0)
            status = stdscr.derwin(1, w, h-1, 0)

            # Dibujar todo de forma segura
            safe_addstr(header, 0, 0, "="*w, curses.A_REVERSE | curses.color_pair(1))
            header.refresh()

            menu.erase()
            safe_addstr(menu, 0, 2, "LAUNDRY DEMONOLOGY TERMINAL v10", curses.A_BOLD|curses.color_pair(1))
            case = random.choices(["YELLOW/AMBER","RED/VERMILLION"],[94,6])[0]
            safe_addstr(menu, 2, 2, f"CASE NIGHTMARE GREEN: {case}",
                       curses.color_pair(3)|curses.A_BOLD if "RED" in case else curses.color_pair(4))
            safe_addstr(menu, 4, 2, "Select operation:", curses.A_BOLD)
            for i, (k, (name, algo, _)) in enumerate(OPTIONS.items(), 6):
                col = curses.color_pair(3)|curses.A_BOLD if k=="4" else curses.color_pair(2)
                safe_addstr(menu, i, 4, f"[{k}] {name:<12} -> {algo}", col)
            safe_addstr(menu, 12, 4, "[Q] Logout & destroy evidence", curses.color_pair(1))
            menu.refresh()

            safe_addstr(sep1, 0, 0, "-"*w, curses.A_DIM); sep1.refresh()
            safe_addstr(sep2, 0, 0, "-"*w, curses.A_DIM); sep2.refresh()
            safe_addstr(warn, 0, 0, random.choice(WARNINGS), curses.color_pair(4)|curses.A_BOLD); warn.refresh()
            safe_addstr(status, 0, 0,
                f"User: {current_user} | Level: {current_clearance} | Entropy: {random.randint(42,97)}%")
            status.refresh()

            if threading.active_count() < 6:
                threading.Thread(target=draw_symbols, args=(logwin,), daemon=True).start()
            print_log(logwin, log_lines)

            key = stdscr.getch()
            if key == -1: continue
            try:
                ch = chr(key).upper()
            except:
                continue

            if ch == "Q":
                stop_symbols.set()
                stdscr.clear()
                msg = "LOGOUT SUCCESSFUL - EVIDENCE DESTROYED"
                safe_addstr(stdscr, h//2, (w-len(msg))//2, msg, curses.A_BOLD|curses.color_pair(3))
                stdscr.refresh()
                time.sleep(3)
                return

            if ch in OPTIONS:
                op = OPTIONS[ch]
                log_lines += ["", f">>> {op[0].upper()} ACTIVE", f"    {op[1]}", ""]
                for phase in range(1,9):
                    log_lines.append(f"    Phase {phase}/8 processing...")
                    print_log(logwin, log_lines)
                    time.sleep(0.4)
                result = random.choice([
                    "SUCCESS: Reality stabilized",
                    "MINOR BREACH: Soul contested",
                    "ERROR: The Sleeper stirs",
                    "Timeline altered",
                    "Entity repelled",
                    "Physics now optional",
                ])
                log_lines += ["", result]
                print_log(logwin, log_lines)
                time.sleep(3)

        except:
            time.sleep(0.1)
            continue

if __name__ == "__main__":
    try:
        curses.wrapper(main)
    except:
        print("\nTerminal destroyed. No evidence remains.")
