#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Simulacion_de_Entorno_Hacker.py
#
#  Copyright 2026 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# ---------------------------
# Importación de los módulos
# ---------------------------
import tkinter as tk
from tkinter import scrolledtext
import time
import random

# --- CONFIGURACIÓN DE ESTILO ---
BACKGROUND_COLOR = '#000000'
FOREGROUND_COLOR = '#00ff41'  # Verde brillante (estilo Matrix)
FONT_FAMILY = 'Consolas'      # Fuente monoespaciada esencial
FONT_SIZE = 10
SPEED_MS = 5                  # Velocidad base de actualización (en milisegundos)

# --- SIMULACIÓN DE CONTENIDO DE GENACT ---

def generate_ansible_output(lines=25):
    """Simula la ejecución de tareas Ansible."""
    output = []
    hosts = ['web-server-01', 'db-server-02', 'proxy-edge-03', 'mgmt-node-10']
    tasks = ['Gathering Facts', 'Applying Security Patches', 'Configuring Firewall', 'Restarting Services']
    states = ['ok', 'changed', 'unreachable', 'failed']
    for i in range(lines):
        host = random.choice(hosts)
        task = random.choice(tasks)
        state = random.choice(states)
        color = 'green' if state in ['ok', 'changed'] else 'red'

        line = f"[{time.strftime('%H:%M:%S')}] {host:<15} : "

        if state == 'ok':
            line += f"\033[92mok\033[0m={random.randint(1, 5)} changed=0 unreachable=0 failed=0 tasks={task}"
        elif state == 'changed':
            line += f"\033[93mchanged\033[0m={random.randint(1, 3)} ok=0 unreachable=0 failed=0 tasks={task}"
        elif state == 'unreachable':
            line += f"\033[91munreachable\033[0m=1 msg='Connection refused.' tasks={task}"
        else: # failed
            line += f"\033[91mfailed\033[0m=1 msg='Permission denied.' tasks={task}"

        output.append(line)
    return output

def generate_botnet_output(lines=25):
    """Simula la actividad de una botnet escaneando y atacando IPs."""
    output = []
    prefixes = ['192.168.', '10.0.', '172.16.']
    ports = [22, 80, 443, 3389, 21, 8080]
    actions = ['SCANNING', 'ATTEMPTING', 'PAYLOAD_DROP', 'CONNECTED']
    for i in range(lines):
        ip = prefixes[random.randint(0, 2)] + str(random.randint(1, 254)) + '.' + str(random.randint(1, 254))
        port = random.choice(ports)
        action = random.choice(actions)

        if action == 'SCANNING':
            line = f"[{time.strftime('%H:%M:%S')}] \033[94m{action:<12}\033[0m: Target IP {ip}:{port}"
        elif action == 'ATTEMPTING':
            line = f"[{time.strftime('%H:%M:%S')}] \033[93m{action:<12}\033[0m: Sending exploit to {ip}. Payload: 'zero-day-v3'"
        elif action == 'CONNECTED':
            line = f"[{time.strftime('%H:%M:%S')}] \033[92m{action:<12}\033[0m: New bot added: {ip}. System: Linux Kernel {random.uniform(5.0, 6.5):.2f}"
        else:
            line = f"[{time.strftime('%H:%M:%S')}] \033[91m{action:<12}\033[0m: Dropping payload on {ip}:{port}"

        output.append(line)
    return output

def generate_bruteforce_output(lines=25):
    """Simula un ataque de fuerza bruta."""
    output = []
    target = "SHA256: 7d251d2f..."
    for i in range(lines):
        guess = ''.join(random.choices('abcdefghijklmnopqrstuvwxyz0123456789', k=random.randint(4, 8)))
        attempts = i * random.randint(100, 500)

        if i % 5 == 0 and i > 0:
            line = f"[{time.strftime('%H:%M:%S')}] \033[93m[MATCH_FAIL]\033[0m HASH: {target} Guess: '{guess}' Attempts: {attempts:,}"
        else:
            line = f"[{time.strftime('%H:%M:%S')}] [TRYING] Guess: '{guess}' Attempts: {attempts:,}"

        output.append(line)
    return output

def generate_julia_output(lines=25):
    """Simula una rutina de cálculo de alta velocidad o simulación científica (Julia)."""
    output = []
    tasks = ['Eigenvalue_Computation', 'MonteCarlo_Simulation', 'Torus_Projection', 'DFT_Analysis']
    for i in range(lines):
        task = random.choice(tasks)
        progress = random.randint(1, 99)
        time_elapsed = random.uniform(0.1, 5.0)

        if progress < 90:
            line = f"[{time.strftime('%H:%M:%S')}] \033[96m[JULIA]\033[0m Task: {task:<25} Progress: {progress:02d}% Time: {time_elapsed:.2f}s"
        else:
            line = f"[{time.strftime('%H:%M:%S')}] \033[92m[COMPLETED]\033[0m Task: {task:<25} Result: {random.uniform(1000, 9999):.4f}"

        output.append(line)
    return output

def generate_memdump_output(lines=25):
    """Simula un volcado de memoria (Memory Dump)."""
    output = []
    for i in range(lines):
        address = f"{i * 0x1000:08X}"
        hex_data = ' '.join([f"{random.randint(0, 255):02X}" for _ in range(16)])
        ascii_data = ''.join([random.choice('abcdefghijklmnopqrstuvwxyz0123456789 .!$') if random.random() > 0.1 else '.' for _ in range(16)])

        line = f"0x{address}: {hex_data} | {ascii_data}"
        output.append(line)
    return output

def generate_docker_build_output(lines=25):
    """Simula la construcción de una imagen Docker."""
    output = []
    layers = ['FROM alpine:latest', 'RUN apk add python3', 'COPY requirements.txt', 'RUN pip install -r', 'CMD ["/app/run.sh"]']
    for i in range(lines):
        layer_id = f"{random.randint(100, 999):03d}"
        action = random.choice(['---> Running in', '---> Using cache', '---> Building layer'])

        if 'cache' in action:
            line = f"Step {i+1}/{len(layers)} : \033[94m{layers[i % len(layers)]}\033[0m"
            line += f"\n {layer_id} {action}..."
        else:
            size_mb = random.uniform(0.1, 10.0)
            line = f"Step {i+1}/{len(layers)} : \033[93m{layers[i % len(layers)]}\033[0m"
            line += f"\n {layer_id} {action} Size: {size_mb:.2f}MB"

        output.append(line)
    return output

# Mapeo de modos a funciones generadoras y títulos
MODES = {
    'ansible': {
        'title': '~/bin/genact -m ansible',
        'generator': generate_ansible_output
    },
    'botnet': {
        'title': '~/bin/genact -m botnet',
        'generator': generate_botnet_output
    },
    'bruteforce': {
        'title': '~/bin/genact -m bruteforce',
        'generator': generate_bruteforce_output
    },
    'julia': {
        'title': '~/bin/genact -m julia',
        'generator': generate_julia_output
    },
    'memdump': {
        'title': '~/bin/genact -m memdump',
        'generator': generate_memdump_output
    },
    'docker_build': {
        'title': '~/bin/genact -m docker_build',
        'generator': generate_docker_build_output
    }
}

# Almacena el estado de cada terminal: (Text Widget, Lista de Líneas, Índice Actual)
terminal_states = {}

def apply_color_tags(text_widget, text_line):
    """Aplica colores a palabras clave si se usa formato ANSI (simulado)."""

    # Reemplazar códigos ANSI (simulados) por tags de Tkinter

    # Verde (ok, connected) - \033[92m
    text_line = text_line.replace('\033[92m', '<green_tag>').replace('\033[0m', '</green_tag>')
    # Amarillo (changed, attempting) - \033[93m
    text_line = text_line.replace('\033[93m', '<yellow_tag>').replace('\033[0m', '</yellow_tag>')
    # Azul/Cyan (scanning, layers) - \033[94m / \033[96m
    text_line = text_line.replace('\033[94m', '<blue_tag>').replace('\033[0m', '</blue_tag>')
    text_line = text_line.replace('\033[96m', '<cyan_tag>').replace('\033[0m', '</cyan_tag>')
    # Rojo (failed, unreachable) - \033[91m
    text_line = text_line.replace('\033[91m', '<red_tag>').replace('\033[0m', '</red_tag>')

    # 1. Dividir la línea por los delimitadores de tag
    parts = []
    current_part = ""
    in_tag = False
    current_tag_name = ""

    # Usar un enfoque basado en cadenas para manejar tags
    segments = []
    # Usamos un marcador de fin común para simplificar
    text_line = text_line.replace('</green_tag>', '||END||').replace('</yellow_tag>', '||END||')
    text_line = text_line.replace('</blue_tag>', '||END||').replace('</cyan_tag>', '||END||')
    text_line = text_line.replace('</red_tag>', '||END||')

    # Tags de inicio
    text_line = text_line.replace('<green_tag>', '||START:green||')
    text_line = text_line.replace('<yellow_tag>', '||START:yellow||')
    text_line = text_line.replace('<blue_tag>', '||START:blue||')
    text_line = text_line.replace('<cyan_tag>', '||START:cyan||')
    text_line = text_line.replace('<red_tag>', '||START:red||')

    # Parsear y aplicar
    parts = text_line.split('||')

    # 0: Texto normal, 1: START:color, 2: Texto coloreado, 3: END, 4: Texto normal...

    for i, part in enumerate(parts):
        if not part:
            continue

        if part.startswith('START:'):
            # Es un marcador de inicio de color. El siguiente segmento es el texto coloreado.
            color_name = part.split(':')[1]
            segments.append(('', color_name))
        elif part == 'END':
            # Es un marcador de fin de color.
            pass
        elif i > 0 and parts[i-1].startswith('START:'):
            # Este es el texto que debe ir coloreado (si el anterior fue START:color)
            color_name = parts[i-1].split(':')[1]
            text_widget.insert(tk.END, part, color_name)
        else:
            # Es texto normal (o el inicio de la línea)
            text_widget.insert(tk.END, part)

    text_widget.insert(tk.END, '\n')
    text_widget.see(tk.END) # Asegura el scroll automático


def update_terminal(terminal_id):
    """Función que maneja la animación de escritura en cada terminal."""
    try:
        text_widget, lines, current_index = terminal_states[terminal_id]
    except KeyError:
        return # Si se borró el estado, parar la recursión

    # Si se llegó al final de las líneas, reiniciar la simulación
    if current_index >= len(lines):
        # Limpiar y regenerar contenido
        text_widget.delete('1.0', tk.END)
        new_lines = MODES[terminal_id]['generator'](lines=random.randint(20, 30))
        terminal_states[terminal_id] = (text_widget, new_lines, 0)
        current_index = 0
        lines = new_lines

    # Obtener la línea actual y el widget
    line_to_print = lines[current_index]

    # Imprimir la línea con manejo de color simulado
    apply_color_tags(text_widget, line_to_print)

    # Actualizar el índice y programar la próxima actualización
    terminal_states[terminal_id] = (text_widget, lines, current_index + 1)

    # Programar la próxima actualización con una velocidad aleatoria para el efecto "vivo"
    next_speed = random.randint(SPEED_MS, SPEED_MS + 20)
    root.after(next_speed, lambda: update_terminal(terminal_id))


def create_terminal_window(parent, mode_id, title, generator, row, col):
    """Crea y configura un panel de terminal en la cuadrícula."""

    # 1. Contenedor del marco
    frame = tk.Frame(parent, bg=BACKGROUND_COLOR, padx=5, pady=5)
    frame.grid(row=row, column=col, sticky="nsew", padx=3, pady=3)

    # Hacer que las celdas de la cuadrícula se expandan
    parent.grid_columnconfigure(col, weight=1)
    parent.grid_rowconfigure(row, weight=1)

    # 2. Título (simula el comando ejecutado)
    title_label = tk.Label(frame, text=title, bg=BACKGROUND_COLOR, fg=FOREGROUND_COLOR,
                           font=(FONT_FAMILY, FONT_SIZE + 2, 'bold'), anchor='w')
    title_label.pack(fill='x')

    # 3. Text Widget (la propia terminal)
    text_widget = scrolledtext.ScrolledText(
        frame,
        wrap='word',
        bg=BACKGROUND_COLOR,
        fg=FOREGROUND_COLOR,
        font=(FONT_FAMILY, FONT_SIZE),
        borderwidth=0,
        highlightthickness=0,
        insertbackground=FOREGROUND_COLOR, # Color del cursor
        relief='flat'
    )
    text_widget.pack(fill='both', expand=True)

    # 4. Configurar tags de color (simulación ANSI)
    text_widget.tag_config('green', foreground='#00ff41')
    text_widget.tag_config('yellow', foreground='#ffff00')
    text_widget.tag_config('blue', foreground='#00bfff')
    text_widget.tag_config('cyan', foreground='#00ffff')
    text_widget.tag_config('red', foreground='#ff0000')

    # 5. Inicializar el estado y empezar la simulación
    initial_lines = generator(lines=random.randint(20, 30))
    terminal_states[mode_id] = (text_widget, initial_lines, 0)

    # Iniciar la animación
    update_terminal(mode_id)


# --- INICIALIZACIÓN DE TKINTER ---

# Crear la ventana principal
root = tk.Tk()
root.title("Ambiente Hacker 💻 GENACT Simulation")
root.geometry("1200x800")
root.configure(bg=BACKGROUND_COLOR)

# Distribuir los 6 modos en una cuadrícula de 2x3
modes_list = list(MODES.keys())
row, col = 0, 0

for i, mode_id in enumerate(modes_list):
    mode_data = MODES[mode_id]

    # Asegurarse de que el índice se mantenga dentro de la cuadrícula 2x3
    row = i // 3
    col = i % 3

    create_terminal_window(
        root,
        mode_id,
        mode_data['title'],
        mode_data['generator'],
        row,
        col
    )

# Iniciar el bucle principal de Tkinter
root.mainloop()
