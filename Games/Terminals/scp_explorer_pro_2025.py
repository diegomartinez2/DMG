#!/usr/bin/env python3
# scp_explorer_pro_2025.py
# Fundación SCP – Clasified Archive Explorer PRO – Renderizado Markdown Rico

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
from tkinter import font as tkfont
import os
from pathlib import Path
from PIL import Image, ImageTk
import markdown
from bs4 import BeautifulSoup
import re
import threading
import time
import platform

# =============================================================================
# CONFIGURACIÓN
# =============================================================================
BASE_DIR = Path(".")
VALID_IMAGE_EXTS = (".jpg", ".jpeg", ".png", ".webp", ".gif", ".bmp")
MD_FILE = "data.md"
IMAGE_FILE = "image"

# =============================================================================
# VENTANA PRINCIPAL
# =============================================================================
root = tk.Tk()
root.title("SCP Foundation – Classified Archive Explorer PRO v9.9")
root.configure(bg="#0d1117")

# Maximizar ventana
system = platform.system()
if system == "Windows":
    root.state('zoomed')
else:
    root.attributes('-zoomed', True)

# =============================================================================
# ESTILO
# =============================================================================
style = ttk.Style()
style.theme_use("clam")
style.configure(".", background="#0d1117", foreground="#c9d1d9", font=("Consolas", 11))

# =============================================================================
# HEADER
# =============================================================================
header = tk.Frame(root, bg="#161b22", height=100)
header.pack(fill="x", padx=12, pady=12)
header.pack_propagate(False)

tk.Label(header, text="SCP", font=("Consolas", 52, "bold"), fg="#58a6ff", bg="#161b22").pack(side="left", padx=20)
tk.Label(header, text="CLASSIFIED ARCHIVE · MEMETIC HAZARD DETECTED", font=("Consolas", 18), fg="#8b949e", bg="#161b22").pack(side="left", pady=30)

status_frame = tk.Frame(header, bg="#161b22")
status_frame.pack(side="right", padx=30)
tk.Label(status_frame, text="Site-19 · O5 Council Access", font=("Consolas", 16), fg="#8b949e", bg="#161b22").pack()
alert_label = tk.Label(status_frame, text="ALERT LEVEL: EUCLID", font=("Consolas", 28, "bold"), fg="#00ff88", bg="#161b22")
alert_label.pack()

# =============================================================================
# PANEL IZQUIERDO – Dossiers
# =============================================================================
left = tk.Frame(root, bg="#161b22", width=380)
left.pack(side="left", fill="y", padx=12, pady=(0,12))
left.pack_propagate(False)

tk.Label(left, text="Anomalous Dossiers", font=("Consolas", 15, "bold"), fg="#58a6ff", bg="#161b22").pack(anchor="w", padx=20, pady=(20,10))

listbox = tk.Listbox(left, bg="#0d1117", fg="#c9d1d9", font=("Consolas", 13),
                     selectbackground="#30363d", highlightthickness=0, bd=0, activestyle="none")
listbox.pack(fill="both", expand=True, padx=20, pady=10)

# =============================================================================
# ÁREA CENTRAL – Markdown Rico
# =============================================================================
main_area = tk.Frame(root)
main_area.pack(expand=True, fill="both", padx=12, pady=(0,12))

md_frame = tk.LabelFrame(main_area, text=" Containment Document · Level 4 Clearance Required ",
                         fg="#58a6ff", bg="#161b22", font=("Consolas", 13, "bold"))
md_frame.pack(side="left", fill="both", expand=True, padx=(0,6))

# Text widget con soporte avanzado
md_text = tk.Text(md_frame, bg="#0d1117", fg="#c9d1d9", font=("Segoe UI", 11),
                  wrap="word", relief="flat", bd=0, padx=20, pady=20, spacing1=4, spacing3=8)
md_text.pack(fill="both", expand=True, padx=12, pady=12)

# Scrollbar personalizado
scrollbar = ttk.Scrollbar(md_frame, orient="vertical", command=md_text.yview)
scrollbar.pack(side="right", fill="y")
md_text.config(yscrollcommand=scrollbar.set)

# =============================================================================
# PANEL DERECHO – Imagen principal
# =============================================================================
img_frame = tk.LabelFrame(main_area, text=" Visual Evidence · Item #: IMG-01 ",
                          fg="#58a6ff", bg="#161b22", font=("Consolas", 13, "bold"))
img_frame.pack(side="right", fill="both", expand=True, padx=(6,0))

img_label = tk.Label(img_frame, bg="#0d1117", text="NO VISUAL DATA", fg="#666666", font=("Consolas", 14))
img_label.pack(fill="both", expand=True, padx=12, pady=12)

# =============================================================================
# LOG
# =============================================================================
log_frame = tk.LabelFrame(root, text=" System Log · Cognitive Resistance Active ",
                          fg="#58a6ff", bg="#161b22", font=("Consolas", 13, "bold"))
log_frame.pack(fill="x", padx=12, pady=(0,12))
log_text = scrolledtext.ScrolledText(log_frame, bg="#0d1117", fg="#c9d1d9", font=("Consolas", 10), height=6, state="disabled")
log_text.pack(fill="x", padx=12, pady=8)

# =============================================================================
# RENDERIZADO MARKDOWN RICO
# =============================================================================
def log(msg, level="INFO"):
    log_text.config(state="normal")
    colors = {"INFO":"#c9d1d9", "WARNING":"#ffb02e", "CRITICAL":"#ff4444", "SUCCESS":"#00ff44", "ACCESS":"#58a6ff"}
    log_text.insert(tk.END, f"[{time.strftime('%H:%M:%S')}] {msg}\n", level)
    log_text.tag_config(level, foreground=colors.get(level, "#c9d1d9"))
    log_text.see(tk.END)
    log_text.config(state="disabled")

# Cache de imágenes
image_cache = {}
current_photo = None

def load_image_in_text(img_path, max_width=500):
    if img_path in image_cache:
        return image_cache[img_path]

    try:
        img = Image.open(img_path)
        ratio = max_width / img.width
        new_height = int(img.height * ratio)
        img = img.resize((max_width, new_height), Image.Resampling.LANCZOS)
        photo = ImageTk.PhotoImage(img)
        image_cache[img_path] = photo
        return photo
    except:
        return None

def render_markdown_rich(content, folder_path):
    md_text.config(state="normal")
    md_text.delete(1.0, tk.END)
    md_text.image_create(tk.END, image=None)  # limpiar imágenes anteriores

    if not content.strip():
        md_text.insert(tk.END, "« DOCUMENTO NO ENCONTRADO O CORRUPTO »", "critical")
        md_text.tag_config("critical", foreground="#ff4444", font=("Consolas", 14, "bold"), justify="center")
        md_text.config(state="disabled")
        return

    # Convertir Markdown → HTML
    html = markdown.markdown(content, extensions=['fenced_code', 'tables', 'toc', 'nl2br'])

    # Parsear con BeautifulSoup
    soup = BeautifulSoup(html, 'html.parser')

    for elem in soup.children:
        if elem.name == "h1":
            md_text.insert(tk.END, elem.text + "\n\n", "h1")
        elif elem.name == "h2":
            md_text.insert(tk.END, elem.text + "\n\n", "h2")
        elif elem.name == "h3":
            md_text.insert(tk.END, elem.text + "\n", "h3")
        elif elem.name == "hr":
            md_text.insert(tk.END, "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n", "hr")
        elif elem.name == "blockquote":
            for line in elem.text.split("\n"):
                md_text.insert(tk.END, f"│ {line}\n", "quote")
            md_text.insert(tk.END, "\n")
        elif elem.name == "ul":
            for li in elem.find_all("li", recursive=False):
                md_text.insert(tk.END, f" • {li.get_text()}\n", "bullet")
            md_text.insert(tk.END, "\n")
        elif elem.name == "ol":
            for i, li in enumerate(elem.find_all("li", recursive=False), 1):
                md_text.insert(tk.END, f" {i}. {li.get_text()}\n", "bullet")
            md_text.insert(tk.END, "\n")
        elif elem.name == "pre":
            code_block = elem.find("code")
            lang = code_block.get("class", [""])[0].replace("language-", "") if code_block else ""
            code = code_block.text if code_block else elem.text
            md_text.insert(tk.END, f" {code.strip()}\n\n", "code")
        elif elem.name == "table":
            # Tabla simple (mejorado)
            rows = elem.find_all("tr")
            for row in rows:
                cells = row.find_all(["td", "th"])
                text = " | ".join(cell.get_text(strip=True) for cell in cells)
                tag = "table_header" if row.find("th") else "table_row"
                md_text.insert(tk.END, f" {text}\n", tag)
            md_text.insert(tk.END, "\n")
        elif elem.name == "img":
            src = elem.get("src")
            alt = elem.get("alt", "Imagen")
            img_path = (folder_path / src).resolve()
            if img_path.exists():
                photo = load_image_in_text(img_path)
                if photo:
                    md_text.image_create(tk.END, image=photo, padx=20, pady=10)
                    md_text.insert(tk.END, f"\n[{alt}]\n\n")
            else:
                md_text.insert(tk.END, f"[Imagen no encontrada: {src}]\n\n", "warning")
        elif elem.name == "p":
            text = elem.get_text()
            # Resaltar **negrita**, *cursiva*, __subrayado__
            parts = re.split(r'(\*\*.*?\*\*|\*.*?\*|\_\_.*?\_+|\~\~.*?~~)', text)
            for part in parts:
                if part.startswith("**") and part.endswith("**"):
                    md_text.insert(tk.END, part[2:-2], "bold")
                elif part.startswith("~~") and part.endswith("~~"):
                    md_text.insert(tk.END, part[2:-2], "strikethrough")
                elif part.startswith("*") and part.endswith("*"):
                    md_text.insert(tk.END, part[1:-1], "italic")
                else:
                    md_text.insert(tk.END, part)
            md_text.insert(tk.END, "\n\n")

    # Estilos finales
    md_text.tag_config("h1", font=("Segoe UI", 22, "bold"), foreground="#58a6ff", spacing1=20, spacing3=10)
    md_text.tag_config("h2", font=("Segoe UI", 18, "bold"), foreground="#8b949e", spacing1=15, spacing3=8)
    md_text.tag_config("h3", font=("Segoe UI", 15, "bold"), foreground="#c9d1d9")
    md_text.tag_config("bold", font=("Segoe UI", 11, "bold"), foreground="#ffffff")
    md_text.tag_config("italic", font=("Segoe UI", 11, "italic"), foreground="#88ffff")
    md_text.tag_config("strikethrough", font=("Segoe UI", 11), foreground="#888888", overstrike=True)
    md_text.tag_config("code", background="#161b22", foreground="#79c0ff", font=("Consolas", 10), relief="flat", lmargin1=20, lmargin2=20, spacing1=5, spacing3=5)
    md_text.tag_config("quote", background="#1e1e2e", foreground="#bb9eff", lmargin1=30, lmargin2=30)
    md_text.tag_config("bullet", foreground="#ff8800")
    md_text.tag_config("hr", foreground="#30363d")
    md_text.tag_config("table_header", background="#161b22", foreground="#58a6ff", font=("Consolas", 11, "bold"))
    md_text.tag_config("table_row", background="#0d1117", foreground="#c9d1d9")
    md_text.tag_config("warning", foreground="#ff4444")

    md_text.config(state="disabled")

# =============================================================================
# CARGAR IMAGEN PRINCIPAL
# =============================================================================
def load_main_image(folder_path):
    global current_photo
    img_label.config(text="CLASSIFIED VISUAL DATA", image="")
    for ext in VALID_IMAGE_EXTS:
        img_path = folder_path / f"{IMAGE_FILE}{ext}"
        if img_path.exists():
            try:
                img = Image.open(img_path)
                img.thumbnail((800, 800), Image.Resampling.LANCZOS)
                current_photo = ImageTk.PhotoImage(img)
                img_label.config(image=current_photo, text="")
                log(f"Visual evidence loaded: {img_path.name}", "SUCCESS")
                return
            except Exception as e:
                log(f"Image error: {e}", "WARNING")
    img_label.config(text="NO VISUAL DATA AVAILABLE", image="")

# =============================================================================
# CARGAR CARPETA
# =============================================================================
def load_folder(event=None):
    sel = listbox.curselection()
    if not sel: return
    folder_name = listbox.get(sel[0])
    folder_path = BASE_DIR / folder_name
    if not folder_path.is_dir(): return

    log(f"ACCESSING: {folder_name.upper()}", "ACCESS")
    alert_label.config(text="ALERT LEVEL: EUCLID", fg="#00ff88")

    # Markdown
    md_path = folder_path / MD_FILE
    if md_path.exists():
        try:
            content = md_path.read_text(encoding="utf-8")
            render_markdown_rich(content, folder_path)
            log(f"Document rendered: {MD_FILE} ({len(content.splitlines())} lines)", "SUCCESS")
        except Exception as e:
            log(f"Read error: {e}", "CRITICAL")
    else:
        render_markdown_rich("**ERROR 404:** Containment document not found.", folder_path)
        log("data.md not found", "WARNING")

    # Imagen principal
    threading.Thread(target=load_main_image, args=(folder_path,), daemon=True).start()

# =============================================================================
# INICIO
# =============================================================================
def refresh():
    listbox.delete(0, tk.END)
    folders = [p for p in BASE_DIR.iterdir() if p.is_dir() and not p.name.startswith(".")]
    folders.sort(key=lambda x: x.name.lower())
    for f in folders:
        listbox.insert(tk.END, f.name)
    log(f"Found {len(folders)} anomalous dossiers", "INFO")

log("SCP CLASSIFIED ARCHIVE TERMINAL v9.9 · ONLINE", "SUCCESS")
log("Memetic kill agent: ACTIVE", "WARNING")
log("Cognito-hazard filters: ENGAGED", "WARNING")
refresh()

listbox.bind("<<ListboxSelect>>", load_folder)
listbox.bind("<Double-1>", load_folder)

root.mainloop()
