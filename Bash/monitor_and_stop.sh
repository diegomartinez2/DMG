#!/bin/bash

# Tiempo en segundos entre cada comprobación
INTERVALO=10

CURRENT_USER=$(whoami)

echo "Iniciando monitorización para $CURRENT_USER..."
echo "Comprobando presencia de otros usuarios cada $INTERVALO segundos."

while true; do
    # Contar usuarios ÚNICOS conectados activos (excluyendo al usuario actual)
    OTHER_USERS_COUNT=$(who | awk '{print $1}' | sort -u | grep -v -x "$CURRENT_USER" | wc -l)

    # Si se detecta al menos a un usuario diferente
    if [ "$OTHER_USERS_COUNT" -gt 0 ]; then
        echo "──────────────────────────────────────────────────────────"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ¡AVISO! Se han detectado otros usuarios conectados."
        echo "Pausando el cálculo pesado (lmp_mpi)..."
        echo "──────────────────────────────────────────────────────────"

        # Pausa el proceso lmp_mpi
        pkill -STOP lmp_mpi

        # Rompe el bucle y finaliza el script
        exit 1
    fi

    # Esperar antes de la siguiente verificación
    sleep $INTERVALO
done
