#!/bin/bash
# Script: hacker_session.sh
# Descripción: Crea una sesión de tmux con paneles divididos, cada uno ejecutando
# un modo diferente del simulador de hacking 'genact'.

# --- PRE-REQUISITOS ---
# 1. Asegúrate de tener 'tmux' y 'genact' instalados.
#    - Para tmux: sudo apt install tmux (o equivalente)
#    - Para genact: Busca 'genact' en GitHub para instrucciones de instalación.

# Verifica si tmux está instalado
if ! command -v tmux &> /dev/null
then
    echo "ERROR: tmux no está instalado. Por favor, instálalo para usar este script."
    exit 1
fi

# Verifica si genact está instalado
if ! command -v genact &> /dev/null
then
    echo "ADVERTENCIA: genact no está instalado. Usaremos un comando de simulación básico en su lugar."
    GENACT_CMD="while true; do echo \"[$(date +%H:%M:%S)] SIMULACION: /dev/null is busy. Ping $(shuf -i 1-254 -n 1).$(shuf -i 1-254 -n 1).$(shuf -i 1-254 -n 1).$(shuf -i 1-254 -n 1)\" | lolcat; sleep 0.05; done"
else
    # Si genact existe, usa el binario real
    GENACT_CMD="genact"
fi

SESSION_NAME="HackerSim"

echo "Preparando sesión de simulación hacker..."

# 1. Crear una nueva sesión de tmux (con un nombre único)
tmux new-session -d -s $SESSION_NAME

# 2. Renombrar la ventana inicial
tmux rename-window 'HackerGrid'

# --- CONFIGURACIÓN DE PANELES (GRID 2x3) ---

# El primer panel (pane 0) ya está creado. Ejecutar ansible.
tmux send-keys -t $SESSION_NAME:0.0 "$GENACT_CMD -m ansible" C-m

# 3. Dividir verticalmente el panel 0 para crear el panel 1 (Derecha de Ansible)
tmux split-window -h -t $SESSION_NAME:0.0
tmux send-keys -t $SESSION_NAME:0.1 "$GENACT_CMD -m botnet" C-m

# 4. Dividir verticalmente el panel 1 para crear el panel 2 (Derecha de Botnet)
tmux split-window -h -t $SESSION_NAME:0.1
tmux send-keys -t $SESSION_NAME:0.2 "$GENACT_CMD -m bruteforce" C-m

# 5. Dividir horizontalmente el panel 0 (arriba) para crear el panel 3 (Abajo de Ansible)
tmux split-window -v -t $SESSION_NAME:0.0
tmux send-keys -t $SESSION_NAME:0.3 "$GENACT_CMD -m julia" C-m

# 6. Dividir horizontalmente el panel 1 (arriba) para crear el panel 4 (Abajo de Botnet)
tmux split-window -v -t $SESSION_NAME:0.1
tmux send-keys -t $SESSION_NAME:0.4 "$GENACT_CMD -m memdump" C-m

# 7. Dividir horizontalmente el panel 2 (arriba) para crear el panel 5 (Abajo de Bruteforce)
tmux split-window -v -t $SESSION_NAME:0.2
tmux send-keys -t $SESSION_NAME:0.5 "$GENACT_CMD -m docker_build" C-m

# Opcional: Ajustar el estilo (solo si se usa genact y la terminal lo soporta)
tmux select-pane -t $SESSION_NAME:0.0
tmux set-option default-terminal "screen-256color"

echo "La sesión de tmux '$SESSION_NAME' ha sido creada."
echo "Para unirte a la sesión y ver la simulación, ejecuta:"
echo "tmux attach -t $SESSION_NAME"
echo ""
echo "Para salir y detener la simulación (manteniendo la sesión), presiona Ctrl+B y luego D."
echo "Para cerrar la sesión y detener todos los procesos, usa 'tmux kill-session -t $SESSION_NAME'."
