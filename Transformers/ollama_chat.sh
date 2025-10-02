#!/bin/bash

# Script para iniciar un chat interactivo con Ollama en la terminal
# Soporta Chain-of-Thought (CoT) prompting opcional
# Uso: ./ollama_chat.sh [modelo] [--cot]

set -e

# Configuración
MODEL=${1:-auto}  # Usa 'auto' por defecto, o el modelo especificado
COT=false
if [ "$2" = "--cot" ]; then
    COT=true
fi
OLLAMA_HOST="127.0.0.1:11434"
OLLAMA_BINARY=$(command -v ollama)
TEMP_LOG=$(mktemp)
TEMP_PROMPT=$(mktemp)

# Colores para salida
RED=$(tput setaf 1 2>/dev/null || echo "")
RESET=$(tput sgr0 2>/dev/null || echo "")

# Funciones
error() { echo "${RED}ERROR:${RESET} $*" >&2; exit 1; }
status() { echo ">>> $*" >&2; }

# Verifica si ollama está instalado
if [ -z "$OLLAMA_BINARY" ]; then
    error "Ollama no está instalado. Instálalo con 'conda install ollama' o 'curl -fsSL https://ollama.com/install.sh | sh'."
fi

# Verifica si el servidor está corriendo
check_server() {
    curl -s --fail --connect-timeout 2 http://$OLLAMA_HOST >/dev/null
    return $?
}

# Inicia el servidor si no está corriendo
start_server() {
    status "Iniciando el servidor de Ollama..."
    $OLLAMA_BINARY serve >"$TEMP_LOG" 2>&1 &
    sleep 2  # Espera a que el servidor inicie
    if ! check_server; then
        error "No se pudo iniciar el servidor. Revisa los logs en $TEMP_LOG."
    fi
}

# Verifica si el modelo está descargado
check_model() {
    $OLLAMA_BINARY list | grep -q "$1" 2>/dev/null
    return $?
}

# Descarga el modelo si no está presente
pull_model() {
    status "Descargando el modelo $1..."
    if ! $OLLAMA_BINARY pull "$1"; then
        error "No se pudo descargar el modelo $1."
    fi
}

# Verifica la RAM disponible (en GB)
check_ram() {
    free -m | awk '/Mem:/ {print int($2/1024)}'
}

# Selecciona el modelo según la RAM si se usa 'auto'
select_model() {
    if [ "$MODEL" = "auto" ]; then
        RAM=$(check_ram)
        status "RAM disponible: ${RAM} GB"
        if [ "$RAM" -ge 64 ]; then
            MODEL="llama3.1:70b"
            status "Seleccionando $MODEL (suficiente RAM)"
        elif [ "$RAM" -ge 16 ]; then
            MODEL="llama3:8b"
            status "Seleccionando $MODEL (RAM moderada)"
        else
            error "RAM insuficiente (<16 GB). Se recomienda al menos 16 GB para llama3:8b."
        fi
    fi
}

# Aplica CoT al prompt si está habilitado
apply_cot() {
    local user_prompt="$1"
    if [ "$COT" = true ]; then
        echo "Piensa paso a paso y explica tu razonamiento antes de dar la respuesta final: $user_prompt"
    else
        echo "$user_prompt"
    fi
}

# Limpieza al salir
cleanup() {
    if [ -f "$TEMP_LOG" ]; then
        rm -f "$TEMP_LOG"
    fi
    if [ -f "$TEMP_PROMPT" ]; then
        rm -f "$TEMP_PROMPT"
    fi
}
trap cleanup EXIT

# Lógica principal
status "Verificando el servidor de Ollama..."
if ! check_server; then
    start_server
else
    status "Servidor de Ollama ya está corriendo en $OLLAMA_HOST."
fi

# Selecciona el modelo si es 'auto'
select_model

status "Verificando el modelo $MODEL..."
if ! check_model "$MODEL"; then
    pull_model "$MODEL"
else
    status "El modelo $MODEL ya está descargado."
fi

status "Iniciando chat interactivo con $MODEL..."
if [ "$COT" = true ]; then
    status "Modo Chain-of-Thought activado: las respuestas incluirán razonamiento paso a paso."
else
    status "Modo chat normal activado."
fi
status "Escribe tus mensajes y usa '/exit' o Ctrl+D para salir."

# Bucle interactivo para leer prompts y aplicar CoT
while true; do
    echo -n "> "
    read -r user_input
    if [ "$user_input" = "/exit" ] || [ -z "$user_input" ]; then
        break
    fi
    apply_cot "$user_input" > "$TEMP_PROMPT"
    $OLLAMA_BINARY run "$MODEL" --file "$TEMP_PROMPT"
done
