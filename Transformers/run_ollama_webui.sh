#!/bin/bash

# ==============================================================================
# Script: Ollama + Open-WebUI con Autolimpieza
# ==============================================================================

# --- Configuración ---
CONDA_ENV_NAME="ollama-env"
OLLAMA_COMMAND="ollama serve"
OPEN_WEBUI_COMMAND="open-webui serve"

# --- Variables de control de procesos ---
PID_OLLAMA=""
PID_WEBUI=""

# --- Función de Limpieza (Se ejecuta al salir) ---
cleanup() {
    echo -e "\n\n[!] Iniciando limpieza de procesos..."

    if [ -n "$PID_WEBUI" ]; then
        echo "Deteniendo Open-WebUI (PID: $PID_WEBUI)..."
        kill $PID_WEBUI 2>/dev/null
    fi

    if [ -n "$PID_OLLAMA" ]; then
        echo "Deteniendo Ollama (PID: $PID_OLLAMA)..."
        kill $PID_OLLAMA 2>/dev/null
    fi

    echo "[OK] Todos los servicios detenidos. Saliendo."
    exit 0
}

# Registrar el trap para capturar Ctrl+C (SIGINT) y salidas normales (EXIT)
trap cleanup SIGINT SIGTERM

# --- Inicialización de Conda ---
CONDA_BASE_PATH=$(conda info --base)
source "$CONDA_BASE_PATH/etc/profile.d/conda.sh"

# 1. Activar el entorno
echo "Activando entorno: $CONDA_ENV_NAME..."
conda activate "$CONDA_ENV_NAME" || { echo "Error activando entorno"; exit 1; }

# 2. Iniciar Ollama en segundo plano
echo "Iniciando Ollama..."
$OLLAMA_COMMAND > /dev/null 2>&1 &
PID_OLLAMA=$! # Guardamos el ID del proceso
echo "Ollama iniciado con PID: $PID_OLLAMA"

sleep 2 # Breve espera para estabilidad

# 3. Iniciar Open-WebUI en segundo plano
echo "Iniciando Open-WebUI..."
$OPEN_WEBUI_COMMAND > /dev/null 2>&1 &
PID_WEBUI=$! # Guardamos el ID del proceso
echo "Open-WebUI iniciado con PID: $PID_WEBUI"

# --- Interfaz de usuario ---
echo "------------------------------------------------"
echo " SERVICIOS ACTIVOS"
echo " - Ollama: http://localhost:11434"
echo " - WebUI:  http://localhost:8080 (normalmente)"
echo "------------------------------------------------"
echo "Presiona cualquier tecla para CERRAR todo y salir."
echo "------------------------------------------------"

# Esperar entrada del usuario (bloquea el script aquí)
read -n 1 -s -r

# Al terminar el comando 'read', el script llega al final y dispara el trap de salida
cleanup
