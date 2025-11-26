#!/usr/bin/env python3
# laundry-demonology-terminal.py
# Versión OFICIAL con login fijo + logout

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
    "ALERT: Thaumic turbulence in sub-basement 7",
    "Form 27B/6 required for soul transfers",
    "CASE NIGHTMARE GREEN readiness dropping",
    "The Auditors are watching",
    "R.H.R.: coffee is not a food group",
]

stop_symbols = threading.Event()
current_user = None
current_clearance = None

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
    try:
        win.erase()
        h, w = win.getmaxyx()
        start = max(0, len(lines) - h + 1)
        for i, line in enumerate(lines[start:]):
            if i >= h-1: break
            text = line[:w-1].ljust(w-1)
            if any(x in line.upper() for x in ["ERROR","BREACH","SLEEPER","EVACUATE"]):
                color = curses.A_BOLD | curses.color_pair(3)
            elif "SUCCESS" in line.upper():
                color = curses.A_BOLD | curses.color_pair(2)
            else:
                color = curses.A_NORMAL
            win.addstr(i, 0, text, color)
        win.refresh()
    except:
        pass

def login_screen(stdscr):
    global current_user, current_clearance
    curses.curs_set(1)
    stdscr.nodelay(False)
    stdscr.clear()

    banner = [
        "",
        "    ===========================================",
        "    ==   THE LAUNDRY - OCCULT COMPUTING     ==",
        "    ==   COMPUTATIONAL DEMONOLOGY TERMINAL   ==",
        "    ===========================================",
        "",
        "    CLASSIFIED LEVEL 3 AND ABOVE ONLY",
        "",
    ]
    for i, line in enumerate(banner):
        stdscr.addstr(i+2, 5, line, curses.A_BOLD | curses.color_pair(1))

    stdscr.addstr(12, 5, "Username: ", curses.color_pair(4))
    stdscr.addstr(14, 5, "Clearance (2-4): ", curses.color_pair(4))

    curses.echo()
    user = stdscr.getstr(12, 15, 20).decode().strip() or "guest"
    clearance_input = stdscr.getstr(14, 25, 1).decode().strip()
    curses.noecho()
    curses.curs_set(0)

    try:
        clearance = int(clearance_input)
        if clearance not in (2,3,4):
            clearance = 2
    except:
        clearance = 2

    current_user = user.upper()[:12]
    current_clearance = clearance

    stdscr.addstr(17, 5, "Authentication successful. Loading terminal...", curses.color_pair(2))
    stdscr.refresh()
    time.sleep(2)

def logout_screen(stdscr):
    h, w = stdscr.getmaxyx()
    stdscr.clear()
    lines = [
        "",
        "    LOGOUT SUCCESSFUL",
        f"    User {current_user} session terminated",
        "    All evidence destroyed",
        "    Have a pleasantly apocalyptic day",
        "",
        "    (Terminal self-destruct in 5... 4... 3...)",
    ]
    for i, line in enumerate(lines):
        stdscr.addstr(h//2-4+i, (w-len(line))//2, line, curses.A_BOLD | curses.color_pair(3))
    stdscr.refresh()
    time.sleep(4)

def main(stdscr):
    global stop_symbols, current_user, current_clearance

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

    while True:
        try:
            h, w = stdscr.getmaxyx()
            if h < 24 or w < 70:
                stdscr.erase()
                stdscr.addstr(0,0,"TERMINAL TOO SMALL", curses.A_BOLD)
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

            # Header
            header.addstr(0,0,"="*w, curses.A_REVERSE | curses.color_pair(1))
            header.refresh()

            # Menú
            menu.erase()
            menu.addstr(0,2,"LAUNDRY COMPUTATIONAL DEMONOLOGY v10", curses.A_BOLD|curses.color_pair(1))
            case = random.choices(["YELLOW/AMBER","RED/VERMILLION"],[93,7])[0]
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

            warn.addstr(0,0,random.choice(WARNINGS), curses.color_pair(4)|curses.A_BOLD)
            warn.refresh()

            status.addstr(0,0,
                f"User: {current_user} | Clearance: LEVEL {current_clearance} | Entropy: {random.randint(42,97)}%",
                curses.A_DIM)
            status.refresh()

            # Fondo de símbolos en zona log
            if threading.active_count() < 4:
                threading.Thread(target=draw_symbols, args=(logwin,), daemon=True).start()

            print_log(logwin, log_lines)

            stdscr.nodelay(True)
            key = stdscr.getch()
            stdscr.nodelay(True)

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
                log_lines.append("")
                log_lines.append(f">>> EXECUTING {op[0].upper()}")
                log_lines.append(f"    {op[1]}")
                print_log(logwin, log_lines)

                for phase in range(1,9):
                    log_lines.append(f"    Phase {phase}/8: {random.choice(['Necromantic backprop','Transfinite collapse','Hyperfold'])}...")
                    print_log(logwin, log_lines)
                    time.sleep(0.35 + random.random()*0.45)

                result = random.choice([
                    "SUCCESS: Reality stabilized",
                    "MINOR BREACH: Soul ownership contested",
                    "ERROR: The Sleeper stirs",
                    "JKR inversion successful - timeline altered",
                    "Entity repelled - Angleton owes you one",
                    "Planck constant detuned - physics now optional",
                ])
                log_lines.append(""); log_lines.append(result)
                print_log(logwin, log_lines)
                time.sleep(3)

        except:
            continue

if __name__ == "__main__":
    try:
        curses.wrapper(main)
    except KeyboardInterrupt:
        stop_symbols.set()
        print("\nSession killed. Please report to Decontamination.")
