#!/usr/bin/env python3
# laundry-demonology-terminal.py
# VERSIÓN OFICIAL - THE LAUNDRY COMPUTATIONAL DEMONOLOGY TERMINAL
# Con login, panel dividido, efecto FTAGHU apocalíptico + BEEP en Linux/macOS
# ¡NO EJECUTES EL 4 UN VIERNES!
"""
Este programa simula una terminal de "The Laundry", la oficina del MI6 que mantiene
a raya a los entes sobrenaturales y trata de evitar el fin del mundo...
Simula una terminal de un ordenador de los '90.
"""
import curses
import time
import random
import threading
import sys

# =============================================================================
# CONFIGURACIÓN GLOBAL
# =============================================================================
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
    "Residual Human Resources: lunch is a privilege, not a right",
]

# Variables globales compartidas entre hilos
stop_symbols = threading.Event()      # Para detener el hilo de ruido al salir
noise_window = None                   # La ventana derecha donde se dibuja el ruido
ftaghu_active = False                 # ¿Está activo el modo apocalíptico?
ftaghu_end_time = 0                   # Cuándo termina el efecto Ftaghu
current_user = "GUEST"
current_clearance = 2

# =============================================================================
# FUNCIONES DE UTILIDAD
# =============================================================================
def safe_addstr(win, y, x, text, attr=0):
    """Escribe de forma segura evitando errores si la ventana es demasiado pequeña"""
    if win is None: return
    try:
        h, w = win.getmaxyx()
        win.addstr(y, x, str(text)[:w-x-1], attr)
    except:
        pass

def beep():
    """Emite un beep en terminales que lo soporten (Linux/macOS)"""
    print("\a", end="", flush=True)

# =============================================================================
# HILO DEL RUIDO DEMONÍACO (panel derecho)
# =============================================================================
def draw_noise():
    normal_chars = ".-:=+*#%&@$∞λΔΘΨΩ∮∯∰∇∂ℏ∫∬∭∑∏√∛∜"
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

            # ¿Estamos en modo Ftaghu? → rojo, rápido y caótico
            if ftaghu_active and time.time() < ftaghu_end_time:
                chars = ftaghu_chars
                attr = curses.color_pair(3) | curses.A_BOLD
                speed = 0.08
                # Beep cada segundo durante el caos
                if int(time.time()) != int(time.time() - 0.1):
                    threading.Thread(target=beep, daemon=True).start()
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

# =============================================================================
# IMPRESIÓN DEL LOG (panel izquierdo)
# =============================================================================
def print_log(win, lines):
    if win is None: return
    win.erase()
    h, w = win.getmaxyx()
    for i, line in enumerate(lines[-(h-1):]):
        if i >= h-1: break
        text = str(line)[:w-1]
        if "FTAGHU" in text or "SLEEPER" in text or "NYARLATHOTEP" in text:
            attr = curses.A_BOLD | curses.color_pair(3)
        elif any(x in text.upper() for x in ["ERROR","BREACH"]):
            attr = curses.A_BOLD | curses.color_pair(3)
        elif "SUCCESS" in text.upper():
            attr = curses.A_BOLD | curses.color_pair(2)
        else:
            attr = curses.A_NORMAL
        safe_addstr(win, i, 0, text, attr)
    win.refresh()

# =============================================================================
# PANTALLA DE LOGIN
# =============================================================================
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
    stdscr.move(9, 18);  user_input = stdscr.getstr(9, 18, 20).decode(errors='ignore').strip()
    stdscr.move(10, 25); cl_input = stdscr.getstr(10, 25, 1).decode(errors='ignore').strip()
    curses.noecho()
    curses.curs_set(0)

    current_user = (user_input or "guest").upper()[:12]
    current_clearance = int(cl_input) if cl_input.isdigit() and 2 <= int(cl_input) <= 4 else 2

    safe_addstr(stdscr, 13, 5, "Authentication successful. Loading terminal...", curses.color_pair(2))
    stdscr.refresh()
    time.sleep(1.5)

# =============================================================================
# BUCLE PRINCIPAL - DIBUJO DE LA INTERFAZ
# =============================================================================
def main(stdscr):
    global noise_window, stop_symbols, ftaghu_active, ftaghu_end_time

    # Inicializar colores
    curses.start_color()
    curses.init_pair(1, curses.COLOR_CYAN,    curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_GREEN,   curses.COLOR_BLACK)
    curses.init_pair(3, curses.COLOR_RED,     curses.COLOR_BLACK)
    curses.init_pair(4, curses.COLOR_YELLOW,  curses.COLOR_BLACK)

    stdscr.nodelay(True)
    login_screen(stdscr)

    log_lines = [
        "*** SYSTEM ONLINE ***",
        f"User {current_user} authenticated - LEVEL {current_clearance}",
        "Wards active. Awaiting command.",
    ]

    # Iniciar el hilo del ruido demoníaco (solo una vez)
    threading.Thread(target=draw_noise, daemon=True).start()

    while True:
        try:
            h, w = stdscr.getmaxyx()

            # Protección contra terminal muy pequeño
            if h < 26 or w < 80:
                stdscr.erase()
                safe_addstr(stdscr, 0, 0, "TERMINAL TOO SMALL - RESIZE TO AT LEAST 80x26", curses.A_BOLD|curses.color_pair(3))
                stdscr.refresh()
                time.sleep(1)
                continue

            # =================================================================
            # CREACIÓN DE VENTANAS (derwin = ventana hija)
            # derwin(altura, anchura, y_origen, x_origen)
            # =================================================================
            header       = stdscr.derwin(1,   w,      0,  0)   # Línea 0
            menu         = stdscr.derwin(14,  w,      1,  0)   # Líneas 1 a 14
            sep1         = stdscr.derwin(1,   w,     15,  0)   # Separador
            left_log     = stdscr.derwin(h-19,w//2,  16,  0)   # ← Panel izquierdo
            right_noise  = stdscr.derwin(h-19,w-w//2,16, w//2) # ← Panel derecho (ruido)
            sep2         = stdscr.derwin(1,   w,   h-3,  0)
            warn         = stdscr.derwin(1,   w,   h-2,  0)
            status       = stdscr.derwin(1,   w,   h-1,  0)

            # Actualizamos la ventana que usa el hilo del ruido
            noise_window = right_noise

            # =================================================================
            # DIBUJADO DE CADA SECCIÓN
            # =================================================================
            safe_addstr(header, 0, 0, "═"*w, curses.A_REVERSE | curses.color_pair(1))
            header.refresh()

            menu.erase()
            safe_addstr(menu, 0, 2, "LAUNDRY COMPUTATIONAL DEMONOLOGY TERMINAL v10", curses.A_BOLD|curses.color_pair(1))
            case = random.choices(["YELLOW/AMBER", "RED/VERMILLION"], [94, 6])[0]
            safe_addstr(menu, 2, 2, f"CASE NIGHTMARE GREEN: {case}",
                       curses.color_pair(3)|curses.A_BOLD if "RED" in case else curses.color_pair(4))
            safe_addstr(menu, 4, 2, "Select operation:", curses.A_BOLD)
            for i, (k, (name, algo, _)) in enumerate(OPTIONS.items(), 6):
                col = curses.color_pair(3)|curses.A_BOLD if k == "4" else curses.color_pair(2)
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

            # =================================================================
            # ENTRADA DE TECLADO
            # =================================================================
            key = stdscr.getch()
            if key == -1: continue
            try:
                ch = chr(key).upper()
            except:
                continue

            # Salir
            if ch == "Q":
                stop_symbols.set()
                stdscr.clear()
                msg = "LOGOUT SUCCESSFUL - ALL EVIDENCE DESTROYED"
                safe_addstr(stdscr, h//2, (w-len(msg))//2, msg, curses.A_BOLD|curses.color_pair(3))
                stdscr.refresh()
                time.sleep(3)
                return

            # Ejecutar operación
            if ch in OPTIONS:
                op = OPTIONS[ch]
                log_lines += ["", f">>> EXECUTING {op[0].upper()}", f"    {op[1]}", ""]

                if ch == "4":
                    ftaghu_active = True
                    ftaghu_end_time = time.time() + 15
                    log_lines += [
                        "!!! FTAGHU SEQUENCE INITIATED !!!",
                        "THE SLEEPER STIRS...",
                        "REALITY INTEGRITY: CRITICAL",
                        "EVACUATION ADVISED"
                    ]

                for phase in range(1, 9):
#                    log_lines.append(f"    Phase {phase}/8: Necromantic convergence in progress...")
                    log_lines.append(f"    Phase {phase}/8: {op[2]} in progress...")
                    print_log(left_log, log_lines)
                    time.sleep(0.5)

                if ch == "4":
                    result = random.choice([
                        "THE SLEEPER HAS OPENED ONE EYE",
                        "REALITY BREACH - EVACUATE LONDON IMMEDIATELY",
                        "NYARLATHOTEP SMILES UPON YOU",
                        "ERROR 0xDEADBEEF: SOUL NOT FOUND",
                        "WELCOME TO THE NEW WORLD ORDER",
                    ])
                else:
                    result = random.choice([
                        "SUCCESS: Local reality stabilized",
                        "MINOR BREACH: Soul ownership contested",
                        "Timeline altered successfully",
                        "Entity repelled. Angleton owes you a pint",
                        "Physics now optional in this sector",
                    ])
                log_lines += ["", result]
                print_log(left_log, log_lines)
                time.sleep(5 if ch == "4" else 3)

        except:
            time.sleep(0.1)
            continue

# =============================================================================
# INICIO DEL PROGRAMA
# =============================================================================
if __name__ == "__main__":
    try:
        curses.wrapper(main)
    except KeyboardInterrupt:
        stop_symbols.set()
        print("\n\nEmergency shutdown. Please report to Decontamination Room 101.")
        sys.exit(0)
