#!/bin/bash

# Script para activar un entorno conda y ejecutar ollama serve y open-webui serve

# Propósito: Este script activa un entorno conda específico y ejecuta los comandos necesarios para iniciar ollama serve y open-webui serve.
# Aplica KISS: El script es simple, conciso y fácil de entender.
# Aplica SOLID:
# - Single Responsibility Principle: El script tiene una única responsabilidad: activar el entorno conda y ejecutar los comandos.
# - Open/Closed Principle: El script es abierto para extensión (podrías añadir más comandos fácilmente) pero cerrado para modificación directa del flujo principal.

# Variables
CONDA_ENV_NAME="ollama-env" # Nombre del entorno conda
OLLAMA_COMMAND="ollama serve"
OPEN_WEBUI_COMMAND="open-webui serve"

# Función para verificar si un comando se ejecuta correctamente
run_command() {
  echo "Ejecutando: $1"
  $1
  if [ $? -eq 0 ]; then
    echo "$1 se ejecutó correctamente."
  else
    echo "$1 falló. Salida de error: $(history | tail -n 1)"
    exit 1 # Sale del script si un comando falla
  fi
}

# 1. Activar el entorno conda
echo "Activando el entorno conda: $CONDA_ENV_NAME"
run_command "conda activate $CONDA_ENV_NAME"

# 2. Iniciar ollama serve
echo "Iniciando ollama serve..."
run_command "$OLLAMA_COMMAND"

# 3. Iniciar open-webui serve
echo "Iniciando open-webui serve..."
run_command "$OPEN_WEBUI_COMMAND"
