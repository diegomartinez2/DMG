#!/usr/bin/env python3
# scp_explorer_2025.py
# Fundación SCP – Explorador de Archivos Clasificados (Estilo SCP)
# Navega carpetas que contengan: data.md + image.(jpg|png|webp|gif)

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
import os
import markdown
from pathlib import Path
from PIL import Image, ImageTk
import threading
import time
import platform

# =============================================================================
# CONFIGURACIÓN
# =============================================================================
BASE_DIR = Path(".")  # Directorio actual
VALID_IMAGE_EXTS = (".jpg", ".jpeg", ".png", ".webp", ".gif", ".bmp")
MD_FILE = "data.md"
IMAGE_FILE = "image"  # Buscará image.jpg, image.png, etc.

# =============================================================================
# VENTANA PRINCIPAL
# =============================================================================
root = tk.Tk()
root.title("SCP Classified Archive Explorer – Secure Access Terminal v9.1")
root.configure(bg="#0d1117")

# Maximizar ventana según SO
system = platform.system()
if system == "Windows":
    root.state('zoomed')
elif system == "Darwin":
    root.attributes('-zoomed', True)
else:
    root.attributes('-zoomed', True)

# Estilo oscuro
style = ttk.Style()
style.theme_use("clam")
style.configure(".", background="#0d1117", foreground="#c9d1d9", font=("Consolas", 11))
style.configure("TLabel", background="#0d1117", foreground="#c9d1d9")

# =============================================================================
# HEADER
# =============================================================================
header = tk.Frame(root, bg="#161b22", height=100)
header.pack(fill="x", padx=12-pady=12)
header.pack_propagate(False)

tk.Label(header, text="SCP", font=("Consolas", 52, "bold"), fg="#58a6ff", bg="#161b22").pack(side="left", padx=20)
tk.Label(header, text="CLASSIFIED ARCHIVE EXPLORER", font=("Consolas", 18), fg="#8b949e", bg="#161b22").pack(side="left", pady=30)

status_frame = tk.Frame(header, bg="#161b22")
status_frame.pack(side="right", padx=30)
tk.Label(status_frame, text="Site-19 Archive", font=("Consolas", 16), fg="#8b949e", bg="#161b22").pack()
alert_label = tk.Label(status_frame, text="ALERT LEVEL: GREEN", font=("Consolas", 28, "bold"), fg="#00ff00", bg="#161b22")
alert_label.pack()

# =============================================================================
# PANEL IZQUIERDO – Lista de carpetas
# =============================================================================
left = tk.Frame(root, bg="#161b22", width=380)
left.pack(side="left", fill="y", padx=12, pady=(0,12))
left.pack_propagate(False)

tk.Label(left, text="Classified Dossiers", font=("Consolas", 15, "bold"), fg="#58a6ff", bg="#161b22").pack(anchor="w", padx=20, pady=(20,10))

listbox = tk.Listbox(left, bg="#0d1117", fg="#c9d1d9", font=("Consolas", 13),
                     selectbackground="#30363d", highlightthickness=0, bd=0)
listbox.pack(fill="both", expand=True, padx=20, pady=10)

# =============================================================================
# ÁREA CENTRAL Y DERECHA
# =============================================================================
main_area = tk.Frame(root)
main_area.pack(expand=True, fill="both", padx=12, pady=(0,12))

# Panel central: Markdown
md_frame = tk.LabelFrame(main_area, text=" Containment Document (data.md) ", fg="#58a6ff", bg="#161b22", font=("Consolas", 13, "bold"))
md_frame.pack(side="left", fill="both", expand=True, padx=(0,6))

md_text = scrolledtext.ScrolledText(md_frame, bg="#0d1117", fg="#c9d1d9", font=("Consolas", 11), wrap="word")
md_text.pack(fill="both", expand=True, padx=12, pady=12)
md_text.tag_config("h1", font=("Consolas", 20, "bold"), foreground="#58a6ff")
md_text.tag_config("h2", font=("Consolas", 16, "bold"), foreground="#8b949e")
md_text.tag_config("bold", font=("Consolas", 11, "bold"))
md_text.tag_config("italic", font=("Consolas", 11, "italic"))

# Panel derecho: Imagen
img_frame = tk.LabelFrame(main_area, text=" Visual Documentation ", fg="#58a6ff", bg="#161b22", font=("Consolas", 13, "bold"))
img_frame.pack(side="right", fill="both", expand=True, padx=(6,0))

img_label = tk.Label(img_frame, bg="#0d1117", text="No image available", fg="#666")
img_label.pack(fill="both", expand=True, padx=12, pady=12)

# Log inferior
log_frame = tk.LabelFrame(root, text=" System Log ", fg="#58a6ff", bg="#161b22", font=("Consolas", 13, "bold"))
log_frame.pack(fill="x", padx=12, pady=(0,12))
log_text = scrolledtext.ScrolledText(log_frame, bg="#0d1117", fg="#c9d1d9", font=("Consolas", 10), height=6, state="disabled")
log_text.pack(fill="x", padx=12, pady=8)

# =============================================================================
# FUNCIONES
# =============================================================================
def log(msg, level="INFO"):
    log_text.config(state="normal")
    colors = {"INFO":"#c9d1d9", "WARNING":"#ffb02e", "CRITICAL":"#ff4444", "SUCCESS":"#00ff44", "ACCESS":"#58a6ff"}
    timestamp = time.strftime('%H:%M:%S')
    log_text.insert(tk.END, f"[{timestamp}] {msg}\n", level)
    log_text.tag_config(level, foreground=colors.get(level, "#c9d1d9"))
    log_text.see(tk.END)
    log_text.config(state="disabled")

def render_markdown(text):
    md_text.config(state="normal")
    md_text.delete(1.0, tk.END)

    if not text.strip():
        md_text.insert(tk.END, "No data.md found in this directory.")
        md_text.config(state="disabled")
        return

    html = markdown.markdown(text, extensions=['fenced_code', 'tables', 'toc'])

    # Convertir HTML simple a texto con formato básico
    lines = text.splitlines()
    in_code = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("```"):
            in_code = not in_code
            md_text.insert(tk.END, line + "\n", "code" if in_code else "")
            continue
        if stripped.startswith("# "):
            md_text.insert(tk.END, line[2:] + "\n", "h1")
        elif stripped.startswith("## "):
            md_text.insert(tk.END, line[3:] + "\n", "h2")
        elif stripped.startswith("**") or stripped.startswith("__"):
            md_text.insert(tk.END, line + "\n", "bold")
        elif stripped.startswith("* ") or stripped.startswith("- "):
            md_text.insert(tk.END, "  • " + line[2:] + "\n")
        else:
            md_text.insert(tk.END, line + "\n")

    md_text.tag_config("code", background="#161b22", font=("Consolas", 10))
    md_text.config(state="disabled")

def load_image(folder_path):
    global current_photo
    img_label.config(text="Loading image...", image="")

    for ext in VALID_IMAGE_EXTS:
        img_path = folder_path / f"{IMAGE_FILE}{ext}"
        if img_path.exists():
            try:
                img = Image.open(img_path)
                # Redimensionar manteniendo proporción
                img.thumbnail((img_label.winfo_width() or 600, img_label.winfo_height() or 800), Image.Resampling.LANCZOS)
                current_photo = ImageTk.PhotoImage(img)
                img_label.config(image=current_photo, text="")
                log(f"Image loaded: {img_path.name}", "SUCCESS")
                return
            except Exception as e:
                log(f"Error loading image: {e}", "WARNING")

    img_label.config(text="No image found\n(image.jpg/png/webp expected)", image="")

def load_folder(event=None):
    selection = listbox.curselection()
    if not selection:
        return
    folder_name = listbox.get(selection[0])
    folder_path = BASE_DIR / folder_name

    if not folder_path.is_dir():
        return

    log(f"Accessing classified dossier: {folder_name}", "ACCESS")

    # Cargar Markdown
    md_path = folder_path / MD_FILE
    if md_path.exists():
        try:
            with open(md_path, "r", encoding="utf-8") as f:
                content = f.read()
            render_markdown(content)
            log(f"Document loaded: data.md ({len(content.splitlines())} lines)", "SUCCESS")
        except Exception as e:
            log(f"Error reading data.md: {e}", "CRITICAL")
    else:
        render_markdown("")
        log("data.md not found in this directory", "WARNING")

    # Cargar imagen en hilo separado
    threading.Thread(target=load_image, args=(folder_path,), daemon=True).start()

# =============================================================================
# CARGAR CARPETAS AL INICIO
# =============================================================================
def refresh_folders():
    listbox.delete(0, tk.END)
    folders = [p for p in BASE_DIR.iterdir() if p.is_dir() and not p.name.startswith(".")]
    folders.sort(key=lambda x: x.name.lower())

    for folder in folders:
        listbox.insert(tk.END, folder.name)

    log(f"Scanned {len(folders)} classified dossiers in current directory", "INFO")

# =============================================================================
# INICIO
# =============================================================================
current_photo = None

log("Secure Archive Terminal v9.1 – ONLINE", "SUCCESS")
log("Memetic kill agent: ARMED", "WARNING")
log("Scanning for classified dossiers...", "INFO")

refresh_folders()

# Bindings
listbox.bind("<<ListboxSelect>>", load_folder)
listbox.bind("<Double-1>", load_folder)

# Iniciar
root.mainloop()
