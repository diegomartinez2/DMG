#!/usr/bin/env python3
# laundry-demonology-terminal.py
# VERSIÓN FINAL CORREGIDA - Panel derecho visible + Ftaghu rojo

import curses
import time
import random
import threading

OPTIONS = {
    "1": ("Protection", "Rosencraft toroidal sequence", "Slow ectastic entity advance"),
    "2": ("Summon", "Level-7 imaginary prime series", "Attract mesonic-energy entities"),
    "3": ("Cast", "nth-cardinality polynomial set", "Disturb JKR boson fields"),
    "4": ("Ftaghu", "Quintic reciprocal convergence", "WAKE THE SLEEPER – DO NOT"),
    "5": ("Expel", "Golberg entirion matrices", "Repel boundary anomalies"),
}

WARNINGS = [
    "WARNING: Unauthorized use violates Benthic Treaty",
    "ALERT: Thaumic turbulence in sub-basement",
    "Form 27B/6 required for soul transfers",
    "CASE NIGHTMARE GREEN readiness dropping",
    "The Auditors are watching",
    "DO NOT RUN FTAGHU ON FRIDAY",
]

stop_symbols = threading.Event()
noise_window = None           # ← Ventana que usará el hilo
ftaghu_active = False
ftaghu_end_time = 0

def safe_addstr(win, y, x, text, attr=0):
    if win is None: return
    try:
        h, w = win.getmaxyx()
        win.addstr(y, x, str(text)[:w-x-1], attr)
    except:
        pass

def draw_noise():
    normal_chars = ".-:=+*#%&@$∞λΔΘΨΩ∮∯∰∇∂ℏ∫∬∭"
    ftaghu_chars = "ΦΨΩΞϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴϴ"

    while not stop_symbols.is_set():
        win = noise_window
        if win is None:
            time.sleep(0.5)
            continue

        try:
            h, w = win.getmaxyx()
            if h < 3 or w < 10:
                time.sleep(0.5); continue

            if ftaghu_active and time.time() < ftaghu_end_time:
                chars = ftaghu_chars
                attr = curses.color_pair(3) | curses.A_BOLD
                speed = 0.08
            else:
                chars = normal_chars
                attr = curses.A_DIM
                speed = 0.35

            win.erase()
            for y in range(h):
                line = "".join(random.choice(chars) for _ in range(w))
                safe_addstr(win, y, 0, line, attr)
            win.refresh()
        except:
            pass
        time.sleep(speed)

def print_log(win, lines):
    if win is None: return
    win.erase()
    h, w = win.getmaxyx()
    for i, line in enumerate(lines[-(h-1):]):
        if i >= h-1: break
        text = str(line)[:w-1]
        if "FTAGHU" in text or "SLEEPER" in text:
            attr = curses.A_BOLD | curses.color_pair(3)
        elif any(x in text.upper() for x in ["ERROR","BREACH"]):
            attr = curses.A_BOLD | curses.color_pair(3)
        elif "SUCCESS" in text.upper():
            attr = curses.A_BOLD | curses.color_pair(2)
        else:
            attr = curses.A_NORMAL
        safe_addstr(win, i, 0, text, attr)
    win.refresh()

def login_screen(stdscr):
    global current_user, current_clearance
    stdscr.clear()
    curses.curs_set(1)
    stdscr.nodelay(False)

    banner = [
        "    ==========================================",
        "    ==   THE LAUNDRY - OCCULT COMPUTING     ==",
        "    ==   COMPUTATIONAL DEMONOLOGY TERMINAL  ==",
        "    ==========================================",
        "    CLASSIFIED LEVEL 3 AND ABOVE ONLY",
        "",
        "    Username:                                      ",
        "    Clearance (2-4):                              ",
    ]
    for i, line in enumerate(banner):
        safe_addstr(stdscr, 3+i, 5, line, curses.A_BOLD | curses.color_pair(1))

    curses.echo()
    stdscr.move(9, 18); user_input = stdscr.getstr(9, 18, 20).decode(errors='ignore').strip()
    stdscr.move(10, 25); cl_input = stdscr.getstr(10, 25, 1).decode(errors='ignore').strip()
    curses.noecho()
    curses.curs_set(0)

    current_user = (user_input or "guest").upper()[:12]
    current_clearance = int(cl_input) if cl_input.isdigit() and 2 <= int(cl_input) <= 4 else 2

    safe_addstr(stdscr, 13, 5, "Authentication successful. Loading...", curses.color_pair(2))
    stdscr.refresh()
    time.sleep(1.5)

def main(stdscr):
    global noise_window, stop_symbols, ftaghu_active, ftaghu_end_time

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

    # Iniciar el hilo del ruido UNA SOLA VEZ
    threading.Thread(target=draw_noise, daemon=True).start()

    while True:
        try:
            h, w = stdscr.getmaxyx()
            if h < 26 or w < 80:
                stdscr.erase()
                safe_addstr(stdscr, 0, 0, "TERMINAL TOO SMALL (min 80x26)", curses.A_BOLD|curses.color_pair(3))
                stdscr.refresh()
                time.sleep(1)
                continue

            # === VENTANAS ===
            header      = stdscr.derwin(1, w, 0, 0)
            menu        = stdscr.derwin(14, w, 1, 0)
            sep1        = stdscr.derwin(1, w, 15, 0)
            left_log    = stdscr.derwin(h-19, w//2, 16, 0)
            right_noise = stdscr.derwin(h-19, w-w//2, 16, w//2)   # ← esta es la buena
            sep2        = stdscr.derwin(1, w, h-3, 0)
            warn        = stdscr.derwin(1, w, h-2, 0)
            status      = stdscr.derwin(1, w, h-1, 0)

            # ACTUALIZAR LA VENTANA QUE USA EL HILO
            noise_window = right_noise

            # Dibujar interfaz
            safe_addstr(header, 0, 0, "═"*w, curses.A_REVERSE | curses.color_pair(1))
            header.refresh()

            menu.erase()
            safe_addstr(menu, 0, 2, "LAUNDRY DEMONOLOGY TERMINAL v10", curses.A_BOLD|curses.color_pair(1))
            case = random.choices(["YELLOW/AMBER","RED/VERMILLION"],[94,6])[0]
            safe_addstr(menu, 2, 2, f"CASE NIGHTMARE GREEN: {case}",
                       curses.color_pair(3)|curses.A_BOLD if "RED" in case else curses.color_pair(4))
            safe_addstr(menu, 4, 2, "Select operation:", curses.A_BOLD)
            for i,(k,(name,algo,_)) in enumerate(OPTIONS.items(),6):
                col = curses.color_pair(3)|curses.A_BOLD if k=="4" else curses.color_pair(2)
                safe_addstr(menu, i, 4, f"[{k}] {name:<12} → {algo}", col)
            safe_addstr(menu, 12, 4, "[Q] Logout & destroy evidence", curses.color_pair(1))
            menu.refresh()

            safe_addstr(sep1, 0, 0, "─"*w, curses.A_DIM); sep1.refresh()
            safe_addstr(sep2, 0, 0, "─"*w, curses.A_DIM); sep2.refresh()
            safe_addstr(warn, 0, 0, random.choice(WARNINGS), curses.color_pair(4)|curses.A_BOLD); warn.refresh()
            safe_addstr(status, 0, 0,
                f"User: {current_user} | Level: {current_clearance} | Entropy: {random.randint(42,97)}%")
            status.refresh()

            print_log(left_log, log_lines)

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
                log_lines += ["", f">>> EXECUTING {op[0].upper()}", f"    {op[1]}", ""]

                if ch == "4":
                    ftaghu_active = True
                    ftaghu_end_time = time.time() + 15
                    log_lines += ["WARNING: FTAGHU SEQUENCE INITIATED",
                                  "THE SLEEPER STIRS...",
                                  "REALITY INTEGRITY: CRITICAL"]

                for phase in range(1,9):
                    log_lines.append(f"    Phase {phase}/8: Necromantic convergence...")
                    print_log(left_log, log_lines)
                    time.sleep(0.5)

                if ch == "4":
                    result = random.choice([
                        "THE SLEEPER HAS OPENED ONE EYE",
                        "REALITY BREACH - EVACUATE LONDON",
                        "NYARLATHOTEP SMILES UPON YOU",
                        "ERROR 0xDEADBEEF: SOUL NOT FOUND",
                    ])
                else:
                    result = random.choice([
                        "SUCCESS: Reality stabilized",
                        "MINOR BREACH: Soul contested",
                        "Timeline altered",
                        "Entity repelled",
                        "Physics now optional",
                    ])
                log_lines += ["", result]
                print_log(left_log, log_lines)
                time.sleep(4 if ch != "4" else 6)

        except:
            time.sleep(0.1)
            continue

if __name__ == "__main__":
    curses.wrapper(main)
