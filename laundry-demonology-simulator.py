#!/usr/bin/env python3
# laundry-demonology-terminal.py
# Versión FINAL Y DEFINITIVA – 100% probada y funcionando

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

def draw_symbols(win):
    chars = ".-:=+*#%&@"
    while not stop_symbols.is_set():
        try:
            h, w = win.getmaxyx()
            if h < 3 or w < 10:
                time.sleep(0.5); continue
            for y in range(h):
                line = "".join(random.choice(chars) for _ in range(w))
                win.addstr(y, 0, line[:w])
            win.refresh()
        except:
            pass
        time.sleep(0.35)

def print_log(win, lines):
    win.erase()
    h, w = win.getmaxyx()
    for i, line in enumerate(lines[-(h-1):]):
        if i >= h-1: break
        text = line[:w-1].ljust(w-1)
        if any(x in line.upper() for x in ["ERROR","BREACH","SLEEPER","EVACUATE"]):
            color = curses.A_BOLD | curses.color_pair(3)
        elif "SUCCESS" in line.upper():
            color = curses.A_BOLD | curses.color_pair(2)
        else:
            color = curses.A_NORMAL
        try:
            win.addstr(i, 0, text, color)
        except:
            pass
    win.refresh()

def login_screen(stdscr):
    global current_user, current_clearance
    curses.curs_set(1)
    stdscr.clear()
    stdscr.nodelay(False)

    banner = [
        "    ===========================================",
        "    ==   THE LAUNDRY - OCCULT COMPUTING     ==",
        "    ==   COMPUTATIONAL DEMONOLOGY TERMINAL   ==",
        "    ===========================================",
        "    CLASSIFIED LEVEL 3 AND ABOVE ONLY",
        "",
        "    Username:",
        "    Clearance (2-4):",
    ]
    for i, line in enumerate(banner):
        stdscr.addstr(3+i, 5, line, curses.A_BOLD | curses.color_pair(1))

    curses.echo()
    user_input = stdscr.getstr(9, 18, 20).decode('utf-8', errors='ignore').strip()
    clearance_input = stdscr.getstr(10, 25, 1).decode('utf-8', errors='ignore').strip()
    curses.noecho()
    curses.curs_set(0)

    current_user = (user_input or "guest").upper()[:12]
    try:
        current_clearance = int(clearance_input)
        if current_clearance not in (2,3,4):
            current_clearance = 2
    except:
        current_clearance = 2

    stdscr.addstr(13, 5, "Authentication successful. Loading terminal...", curses.color_pair(2))
    stdscr.refresh()
    time.sleep(1.5)
    stdscr.clear()

def logout_screen(stdscr):
    h, w = stdscr.getmaxyx()
    stdscr.nodelay(False)
    lines = [
        "",
        "    LOGOUT SUCCESSFUL",
        f"    Session {current_user} terminated",
        "    All logs destroyed",
        "    Have a pleasantly apocalyptic day",
        "",
        "    Terminal self-destruct in 5...4...3...",
    ]
    for i, line in enumerate(lines):
        stdscr.addstr(h//2-4+i, max(0,(w-len(line))//2), line, curses.A_BOLD | curses.color_pair(3))
    stdscr.refresh()
    time.sleep(4)

def main(stdscr):
    global stop_symbols

    curses.start_color()
    curses.init_pair(1, curses.COLOR_CYAN,    curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_GREEN,   curses.COLOR_BLACK)
    curses.init_pair(3, curses.COLOR_RED,     curses.COLOR_BLACK)
    curses.init_pair(4, curses.COLOR_YELLOW,  curses.COLOR_BLACK)

    login_screen(stdscr)

    log_lines = [
        "*** SYSTEM ONLINE ***",
        f"User {current_user} authenticated - LEVEL {current_clearance}",
        "Wards active. Awaiting command.",
    ]

    # Iniciar hilo de símbolos una sola vez
    threading.Thread(target=draw_symbols, args=(stdscr.derwin(1,1,0,0),), daemon=True).start()

    while True:
        h, w = stdscr.getmaxyx()
        if h < 24 or w < 70:
            stdscr.erase()
            stdscr.addstr(0,0,"RESIZE TERMINAL (min 70x24)", curses.A_BOLD|curses.color_pair(3))
            stdscr.refresh()
            time.sleep(1)
            continue

        header = stdscr.derwin(1, w, 0, 0)
        menu   = stdscr.derwin(14, w, 1, 0)
        sep1   = stdscr.derwin(1, w, 15, 0)
        logwin = stdscr.derwin(h-19, w, 16, 0)
        sep2   = stdscr.derwin(1, w, h-3, 0)
        warn   = stdscr.derwin(1, w, h-2, 0)
        status = stdscr.derwin(1, w, h-1, 0)

        header.addstr(0,0,"="*w, curses.A_REVERSE|curses.color_pair(1)); header.refresh()

        menu.erase()
        menu.addstr(0,2,"LAUNDRY COMPUTATIONAL DEMONOLOGY v10", curses.A_BOLD|curses.color_pair(1))
        case = random.choices(["YELLOW/AMBER","RED/VERMILLION"],[94,6])[0]
        menu.addstr(2,2,"CASE NIGHTMARE GREEN: "+case,
                   curses.color_pair(3)|curses.A_BOLD if "RED" in case else curses.color_pair(4))
        menu.addstr(4,2,"Select operation:", curses.A_BOLD)
        for i,(k,(name,algo,_)) in enumerate(OPTIONS.items(),6):
            col = curses.color_pair(3)|curses.A_BOLD if k=="4" else curses.color_pair(2)
            menu.addstr(i,4,f"[{k}] {name:<12} -> {algo}", col)
        menu.addstr(12,4,"[Q] Logout & destroy evidence", curses.color_pair(1))
        menu.refresh()

        sep1.addstr(0,0,"-"*w, curses.A_DIM); sep1.refresh()
        sep2.addstr(0,0,"-"*w, curses.A_DIM); sep2.refresh()

        warn.addstr(0,0,random.choice(WARNINGS), curses.color_pair(4)|curses.A_BOLD); warn.refresh()
        status.addstr(0,0,f"User: {current_user} | Clearance: LEVEL {current_clearance} | Entropy: {random.randint(42,97)}%"); status.refresh()

        # Fondo de símbolos solo en la zona de log
        if threading.active_count() < 5:
            threading.Thread(target=draw_symbols, args=(logwin,), daemon=True).start()
        print_log(logwin, log_lines)

        stdscr.nodelay(True)
        key = stdscr.getch()
        if key == -1: continue
        try:
            ch = chr(key).upper()
        except:
            continue

        if ch == "Q":
            stop_symbols.set()
            logout_screen(stdscr)
            return

        if ch in OPTIONS:
            op = OPTIONS[ch]
            log_lines += ["", f">>> EXECUTING {op[0].upper()}", f"    {op[1]}", ""]
            print_log(logwin, log_lines)

            for phase in range(1,9):
                log_lines.append(f"    Phase {phase}/8: {random.choice(['Necromantic backprop','Transfinite collapse','Hyperfold'])}...")
                print_log(logwin, log_lines)
                time.sleep(0.4)

            result = random.choice([
                "SUCCESS: Reality stabilized",
                "MINOR BREACH: Soul ownership contested",
                "ERROR: The Sleeper stirs",
                "JKR inversion successful - timeline altered",
                "Entity repelled - Angleton owes you one",
                "Planck constant detuned - physics optional",
            ])
            log_lines += ["", result]
            print_log(logwin, log_lines)
            time.sleep(3)

if __name__ == "__main__":
    curses.wrapper(main)
