#!/usr/bin/env python3
# scp_foundation_gui_2025.py
# Fundación SCP – Interfaz Segura 2025 – VERSIÓN UNIVERSAL (Linux/macOS/Windows)

import tkinter as tk
from tkinter import ttk, scrolledtext
import random
import time
import threading
import math
import platform

# =============================================================================
# DATOS DE SCPs
# =============================================================================
SCPS = {
    "SCP-096": ("Shy Guy", "Euclid → Keter si se ve la cara", "#FF4444"),
    "SCP-173": ("The Sculpture", "Euclid – no parpadear", "#FFFF88"),
    "SCP-049": ("Plague Doctor", "Cirugía zombie", "#88FF88"),
    "SCP-106": ("The Old Man", "Keter – corrosión dimensional", "#FF88FF"),
    "SCP-914": ("The Clockworks", "Safe – refinamiento", "#88FFFF"),
    "SCP-079": ("Old AI", "Conciencia digital", "#CCCCCC"),
    "SCP-682": ("Hard-to-Destroy Reptile", "Keter – adaptación total", "#FF0000"),
    "SCP-001": ("[ACCESO DENEGADO]", "Apollyon", "#000000"),
}

ALERTS = ["WHITE", "GREEN", "BLUE", "AMBER", "RED", "KETER"]
current_alert = 2
keter_mode = False

# =============================================================================
# VENTANA PRINCIPAL
# =============================================================================
root = tk.Tk()
root.title("SCP Foundation – Secure Containment Interface v4.7")
root.configure(bg="#0d1117")

# MAXIMIZAR VENTANA SEGÚN EL SISTEMA OPERATIVO
system = platform.system()
if system == "Windows":
    root.state('zoomed')
elif system == "Darwin":  # macOS
    root.attributes('-zoomed', True)
else:  # Linux y otros
    root.attributes('-zoomed', True)   # funciona en la mayoría de WM modernos
    # Alternativa si no funciona lo de arriba:
    # root.geometry(f"{root.winfo_screenwidth()}x{root.winfo_screenheight()-60}+0+0")

# =============================================================================
# ESTILO OSCURO
# =============================================================================
style = ttk.Style()
style.theme_use("clam")
style.configure(".", background="#0d1117", foreground="#c9d1d9", font=("Consolas", 11))
style.configure("TLabel", background="#0d1117", foreground="#c9d1d9")

# =============================================================================
# HEADER
# =============================================================================
header = tk.Frame(root, bg="#161b22", height=100)
header.pack(fill="x", padx=12, pady=12)
header.pack_propagate(False)

tk.Label(header, text="SCP", font=("Consolas", 52, "bold"), fg="#58a6ff", bg="#161b22").pack(side="left", padx=20)
tk.Label(header, text="SECURE · CONTAIN · PROTECT", font=("Consolas", 18), fg="#8b949e", bg="#161b22").pack(side="left", pady=30)

status_frame = tk.Frame(header, bg="#161b22")
status_frame.pack(side="right", padx=30)
tk.Label(status_frame, text="Site-19", font=("Consolas", 16), fg="#8b949e", bg="#161b22").pack()
alert_label = tk.Label(status_frame, text=f"ALERT LEVEL: {ALERTS[current_alert]}",
                      font=("Consolas", 28, "bold"), fg="#00ff00", bg="#161b22")
alert_label.pack()

# =============================================================================
# PANEL IZQUIERDO – Lista de SCPs
# =============================================================================
left = tk.Frame(root, bg="#161b22", width=420)
left.pack(side="left", fill="y", padx=12, pady=(0,12))
left.pack_propagate(False)

tk.Label(left, text="Contained Anomalies", font=("Consolas", 15, "bold"), fg="#58a6ff", bg="#161b22").pack(anchor="w", padx=20, pady=(20,10))

listbox = tk.Listbox(left, bg="#0d1117", fg="#c9d1d9", font=("Consolas", 13),
                    selectbackground="#30363d", highlightthickness=0, bd=0)
for scp in SCPS:
    listbox.insert(tk.END, scp)
listbox.pack(fill="both", expand=True, padx=20, pady=10)

# =============================================================================
# ÁREA CENTRAL Y DERECHA
# =============================================================================
main_area = tk.Frame(root)
main_area.pack(expand=True, fill="both", padx=12, pady=(0,12))

# Log
log_frame = tk.LabelFrame(main_area, text=" Event Log ", fg="#58a6ff", bg="#161b22", font=("Consolas", 13, "bold"))
log_frame.pack(side="left", fill="both", expand=True, padx=(0,6))

log_text = scrolledtext.ScrolledText(log_frame, bg="#0d1117", fg="#c9d1d9", font=("Consolas", 11), state="disabled")
log_text.pack(fill="both", expand=True, padx=12, pady=12)

# Monitor de ondas
wave_frame = tk.LabelFrame(main_area, text=" Anomalous Activity ", fg="#58a6ff", bg="#161b22", font=("Consolas", 13, "bold"))
wave_frame.pack(side="right", fill="both", expand=True, padx=(6,0))

canvas = tk.Canvas(wave_frame, bg="#0d1117", highlightthickness=0)
canvas.pack(fill="both", expand=True, padx=12, pady=12)

# =============================================================================
# FUNCIONES
# =============================================================================
def log(msg, level="INFO"):
    log_text.config(state="normal")
    colors = {"INFO":"#c9d1d9", "WARNING":"#ffb02e", "CRITICAL":"#ff4444", "SUCCESS":"#00ff44"}
    tag = level
    log_text.insert(tk.END, f"[{time.strftime('%H:%M:%S')}] {msg}\n", tag)
    log_text.tag_config(tag, foreground=colors.get(level, "#c9d1d9"))
    log_text.see(tk.END)
    log_text.config(state="disabled")

def beep():
    print("\a", end="", flush=True)

wave_offset = 0
def animate_wave():
    global wave_offset, keter_mode
    if canvas.winfo_width() <= 1:
        root.after(200, animate_wave)
        return

    canvas.delete("all")
    w = canvas.winfo_width()
    h = canvas.winfo_height()
    points = []
    freq = 8 if keter_mode else 2
    amp = 50 if keter_mode else 20
    color = "#ff0000" if keter_mode else "#00ff88"

    for x in range(0, w, 4):
        y = h//2 + math.sin((x + wave_offset) * 0.05 * freq) * amp * random.uniform(0.7, 1.3)
        points.extend([x, y])

    canvas.create_line(points, fill=color, width=3, smooth=True)
    wave_offset += 5
    root.after(50, animate_wave)

def keter_protocol():
    global keter_mode
    keter_mode = True
    alert_label.config(text="ALERT LEVEL: KETER", fg="#ff0000")
    threading.Thread(target=lambda: [beep() or time.sleep(0.4) for _ in range(40)], daemon=True).start()

def on_scp_select(event=None):
    sel = listbox.curselection()
    if not sel: return
    scp = listbox.get(sel[0])

    log(f"Accessing {scp} containment systems...", "INFO")

    if scp == "SCP-096":
        root.after(1200, lambda: log("Image SCP-096-1 detected", "CRITICAL"))
        root.after(2000, lambda: log("Subject breach in T-30 seconds", "CRITICAL"))
    elif scp == "SCP-173":
        root.after(1000, lambda: log("Blink detection: OK", "SUCCESS"))
    elif scp == "SCP-682":
        root.after(1500, lambda: log("Adaptation rate: 412% – containment failing", "CRITICAL"))
        alert_label.config(text="ALERT LEVEL: RED", fg="#ff0000")
    elif scp == "SCP-001":
        keter_protocol()
        log("SCP-001 PROTOCOL ACTIVATED", "CRITICAL")
        log("XK-class end-of-world scenario imminent", "CRITICAL")
    else:
        root.after(1000, lambda: log(f"{scp} containment stable", "SUCCESS"))

# =============================================================================
# BINDINGS
# =============================================================================
listbox.bind("<Double-1>", on_scp_select)
root.bind("<Return>", on_scp_select)

# =============================================================================
# INICIO
# =============================================================================
log("Secure Containment Interface v4.7 – ONLINE", "INFO")
log("Site-19 systems nominal", "SUCCESS")
log("Memetic kill agent: ARMED", "WARNING")

# Animación segura
root.after(500, animate_wave)

root.mainloop()
