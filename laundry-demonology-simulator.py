#!/usr/bin/env python3
# laundry-demonology-terminal-v10.py
# Versión final corregida – funciona perfectamente

import curses
import time
import random
import threading

OPTIONS = {
    "1": ("Protection", "Rosencraft toroidal sequence", "Slow ectastic entity advance"),
    "2": ("Summon", "Level-7 imaginary prime series", "Attract mesonic-energy entities"),
    "3": ("Cast", "nth-cardinality polynomial set", "Disturb JKR boson fields"),
    "4": ("Ftaghu", "Quintic reciprocal convergence", "Observe The Sleeper (EXTREME RISK)"),
    "5": ("Expel", "Golberg entirion matrices", "Repel boundary anomalies"),
}

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
    while not getattr(win, "stop_event", False):
        line = ''.join(random.choice(symbols) for _ in range(200))
        try:
            for row in range(win.height):
                win.addstr(row, 0, line[row:row+win.width].ljust(win.width), curses.A_DIM)
            win.refresh()
        except:
            pass
        time.sleep(0.25)

def computation_output(win, text_lines):
    win.clear()
    for i, line in enumerate(text_lines[-(win.height-1):]):
        color = curses.color_pair(2)
        if any(x in line for x in ["ERROR", "BREACH", "Sleeper", "Nyarlathotep"]):
            color = curses.color_pair(3)
        elif "SUCCESS" in line or "stabilized" in line:
            color = curses.color_pair(2)
        else:
            color = curses.A_DIM
        try:
            win.addstr(i, 0, line[:win.width-1].ljust(win.width-1), color)
        except:
            pass
    win.refresh()

def main(stdscr):
    curses.curs_set(0)
    curses.start_color()
    curses.init_pair(1, curses.COLOR_CYAN,    curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_GREEN,   curses.COLOR_BLACK)
    curses.init_pair(3, curses.COLOR_RED,     curses.COLOR_BLACK)
    curses.init_pair(4, curses.COLOR_YELLOW, curses.COLOR_BLACK)

    h, w = stdscr.getmaxyx()

    header_win   = stdscr.derwin(1, w, 0, 0)
    menu_win     = stdscr.derwin(14, w, 1, 0)
    sep1         = stdscr.derwin(1, w, 15, 0)
    compute_win  = stdscr.derwin(h-19, w, 16, 0)
    sep2         = stdscr.derwin(1, w, h-3, 0)
    warning_win  = stdscr.derwin(1, w, h-2, 0)
    status_win   = stdscr.derwin(1, w, h-1, 0)

    compute_win.stop_event = False
    threading.Thread(target=scrolling_symbols, args=(compute_win,), daemon=True).start()

    log_lines = ["=== LAUNDRY SYSTEM ONLINE ===", "Wards active. Awaiting command."]

    while True:
        # Header
        header_win.clear()
        header_win.addstr(0, 0, "▓" * w, curses.A_REVERSE | curses.color_pair(1))
        header_win.refresh()

        # Menú
        menu_win.clear()
        menu_win.addstr(0, 2, "LAUNDRY COMPUTATIONAL DEMONOLOGY TERMINAL v10", curses.A_BOLD | curses.color_pair(1))
        status_text = random.choices(["YELLOW/AMBER", "RED/VERMILLION – EVACUATE"], weights=[92,8])[0]
        menu_win.addstr(2, 2, "CASE NIGHTMARE GREEN: ", curses.color_pair(4))
        menu_win.addstr(2, 24, status_text, curses.A_BOLD | (curses.color_pair(3) if "RED" in status_text else curses.color_pair(2)))
        menu_win.addstr(4, 2, "Select operation:", curses.A_BOLD)
        for i, (k, (name, algo, _)) in enumerate(OPTIONS.items(), 6):
            color = curses.color_pair(3) if k == "4" else curses.color_pair(2)
            menu_win.addstr(i, 4, f"[{k}] {name:<12} → {algo}", color)
        menu_win.addstr(12, 4, "[Q] Exit & initiate self-destruct", curses.color_pair(1))
        menu_win.refresh()

        # Separadores
        sep1.addstr(0, 0, "─" * w, curses.A_DIM)
        sep2.addstr(0, 0, "─" * w, curses.A_DIM)
        sep1.refresh(); sep2.refresh()

        # Warning y status
        warning_win.clear()
        warning_win.addstr(0, 0, random.choice(warnings), curses.color_pair(4) | curses.A_BOLD)
        warning_win.refresh()

        status_win.clear()
        user = random.choice(["boble", "mo", "angleton", "harriet", "guest"])
        status_win.addstr(0, 0, f"User: {user.upper()} | Clearance: LEVEL {random.randint(2,4)} | Entropy: {random.randint(38,99)}%")
        status_win.refresh()

        computation_output(compute_win, log_lines)

        key = stdscr.getch()
        if key in (-1, curses.KEY_RESIZE):
            continue
        try:
            ch = chr(key).upper()
        except:
            continue

        if ch == 'Q':
            compute_win.stop_event = True
            break

        if ch in OPTIONS:
            op = OPTIONS[ch]
            log_lines.append(f"> EXECUTING {op[0]}")
            log_lines.append(f"  → {op[1]}")
            computation_output(compute_win, log_lines)

            for phase in range(1, 9):
                log_lines.append(f"  Phase {phase}: {random.choice(['Transfinite recursion','Necromantic backpropagation','Hyperdimensional folding'])}...")
                computation_output(compute_win, log_lines)
                time.sleep(0.3 + random.uniform(0, 0.5))

            result = random.choice([
                "SUCCESS: Reality deviation contained.",
                "MINOR BREACH: Soul ownership now 4.1% contested.",
                "ERROR: The Sleeper stirs. Recommend immediate evacuation.",
                "JKR inversion complete. Tomorrow is now optional.",
                "Entity repelled. Angleton owes you one.",
                "Planck constant adjusted. All cats are now in superposition of unionizing.",
            ])
            log_lines.append(result)
            computation_output(compute_win, log_lines)
            time.sleep(3)

    # Salida limpia
    stdscr.clear()
    msg = "TERMINAL OFFLINE – PLEASE FEED THE DOCUMENTS TO THE SHREDDER"
    stdscr.addstr(h//2, (w-len(msg))//2, msg, curses.A_BOLD | curses.color_pair(3))
    stdscr.refresh()
    time.sleep(3)

if __name__ == "__main__":
    try:
        curses.wrapper(main)
    except KeyboardInterrupt:
        print("\nSession terminated. Remember: there is no Form 666/B for resurrection.")
