#!/bin/bash

# ==============================================================================
# Script para activar entorno Conda y ejecutar Ollama + Open-WebUI
# ==============================================================================
# 
# Los scripts de bash no heredan las funciones de shell de conda por defecto.
# Es necesario cargar el script 'conda.sh' de la instalación de base.

# --- Configuración ---
CONDA_ENV_NAME="ollama-env"
OLLAMA_COMMAND="ollama serve"
OPEN_WEBUI_COMMAND="open-webui serve"

# --- Inicialización de Conda dentro del Script ---
# Buscamos la ruta de CONDA_EXE y derivamos la ubicación de conda.sh
CONDA_BASE_PATH=$(conda info --base)
source "$CONDA_BASE_PATH/etc/profile.d/conda.sh"

# Función para verificar ejecución
run_status() {
    if [ $? -eq 0 ]; then
        echo "[OK] $1"
    else
        echo "[ERROR] $1 falló."
        exit 1
    fi
}

# 1. Activar el entorno
echo "Activando entorno: $CONDA_ENV_NAME..."
conda activate "$CONDA_ENV_NAME"
run_status "Activación de entorno"

# 2. Iniciar Ollama en segundo plano (Background)
# Usamos '&' porque 'ollama serve' es un proceso persistente que bloquea el hilo.
echo "Iniciando $OLLAMA_COMMAND en segundo plano..."
$OLLAMA_COMMAND > /dev/null 2>&1 &
sleep 2 # Esperamos un momento a que el servicio levante

# 3. Iniciar Open-WebUI
# Este proceso se mantiene en primer plano para ver los logs y mantener el script vivo.
echo "Iniciando $OPEN_WEBUI_COMMAND..."
$OPEN_WEBUI_COMMAND

# Nota: Si cierras el script con Ctrl+C, el proceso de ollama podría quedar
# corriendo en segundo plano dependiendo de la configuración de tu sistema.
