#!/bin/bash

# Archivos temporales para salida y errores
TEMP_OUT=$(mktemp)
TEMP_ERR=$(mktemp)
TEMP_HARM_PROGRESS=$(mktemp)
TEMP_CUBIC_PROGRESS=$(mktemp)

# Función para limpiar archivos temporales
cleanup() {
    rm -f "$TEMP_OUT" "$TEMP_ERR" "$TEMP_HARM_PROGRESS" "$TEMP_CUBIC_PROGRESS"
}

# Capturar señales para limpieza
trap cleanup EXIT

# Obtener dimensiones de la terminal
TERMINAL_HEIGHT=$(tput lines)
TERMINAL_WIDTH=$(tput cols)

# Calcular dimensiones de las ventanas
OUT_HEIGHT=$((TERMINAL_HEIGHT * 90 / 100 - 2))  # 90% para salida estándar
ERR_HEIGHT=4                                    # Un par de líneas para errores
SEP_HEIGHT=1                                    # Separador
OUT_WIDTH=$((TERMINAL_WIDTH - 10))              # Ancho ajustado
ERR_WIDTH=$((TERMINAL_WIDTH - 10))

# Total de iteraciones (180 + 3348)
TOTAL_ITERS=3528

# Función para procesar un archivo (harm o cubic)
process_file() {
    TYPE=$1  # harm o cubic
    NUM=$2   # número de archivo
    TOTAL=$3 # total de iteraciones para este tipo
    PROGRESS_FILE=$4 # archivo para progreso

    # Formatear número según el tipo
    if [ "$TYPE" = "harm" ]; then
        NUM_FORMATTED=$(printf "%03d" $NUM)
    else
        NUM_FORMATTED=$(printf "%04d" $NUM)
    fi

    # Procesar archivo
    echo "Procesando ${TYPE}${NUM_FORMATTED}.lammps..." >> "$TEMP_OUT"
    cp ${TYPE}${NUM_FORMATTED}.lammps tmp.lammps 2>> "$TEMP_ERR"
    ./lmp_mpi_chimes < in.sw >> "$TEMP_OUT" 2>> "$TEMP_ERR"
    cp XFSET XFSET.${TYPE}${NUM_FORMATTED} 2>> "$TEMP_ERR"

    # Actualizar progreso para este tipo
    PERCENT=$((NUM * 100 / TOTAL))
    echo $PERCENT > "$PROGRESS_FILE"
}

export -f process_file
export TEMP_OUT TEMP_ERR

# Función para monitorear salida y errores
monitor_output() {
    dialog --title "Salida Estándar - LAMMPS" --tailboxbg "$TEMP_OUT" $OUT_HEIGHT $OUT_WIDTH \
        --begin 2 2 &
    DIALOG_PID_OUT=$!

    dialog --title "Errores" --tailboxbg "$TEMP_ERR" $ERR_HEIGHT $ERR_WIDTH \
        --begin $((OUT_HEIGHT + SEP_HEIGHT + 2)) 2 &
    DIALOG_PID_ERR=$!
}

# Función para monitorear progreso total
monitor_total_progress() {
    {
        while [ -f "$TEMP_HARM_PROGRESS" ] || [ -f "$TEMP_CUBIC_PROGRESS" ]; do
            HARM_PERCENT=$(cat "$TEMP_HARM_PROGRESS" 2>/dev/null || echo 0)
            CUBIC_PERCENT=$(cat "$TEMP_CUBIC_PROGRESS" 2>/dev/null || echo 0)
            HARM_COUNT=$((HARM_PERCENT * 180 / 100))
            CUBIC_COUNT=$((CUBIC_PERCENT * 3348 / 100))
            TOTAL_COUNT=$((HARM_COUNT + CUBIC_COUNT))
            TOTAL_PERCENT=$((TOTAL_COUNT * 100 / TOTAL_ITERS))
            echo $TOTAL_PERCENT
            sleep 0.5
        done
    } | dialog --title "Progreso Total" --gauge "Progreso global..." 8 $OUT_WIDTH 0
}

# Iniciar monitoreo de salida y errores
monitor_output

# Iniciar barras de progreso para harm y cubic
{
    seq 1 180 | parallel -j 4 process_file harm {} 180 "$TEMP_HARM_PROGRESS" &
    HARM_PID=$!
    seq 1 3348 | parallel -j 4 process_file cubic {} 3348 "$TEMP_CUBIC_PROGRESS" &
    CUBIC_PID=$!

    # Mostrar barras de progreso para harm y cubic
    dialog --title "Progreso Bucle Harm (1-180)" --gauge "Procesando harm..." 8 $OUT_WIDTH 0 \
        < "$TEMP_HARM_PROGRESS" &
    DIALOG_HARM_PID=$!

    dialog --title "Progreso Bucle Cubic (1-3348)" --gauge "Procesando cubic..." 8 $OUT_WIDTH 0 \
        < "$TEMP_CUBIC_PROGRESS" &
    DIALOG_CUBIC_PID=$!

    # Monitorear progreso total
    monitor_total_progress &

    # Esperar a que terminen los procesos
    wait $HARM_PID $CUBIC_PID
} &

# Esperar a que terminen las barras de progreso
wait $DIALOG_HARM_PID $DIALOG_CUBIC_PID

# Terminar procesos de dialog de salida y errores
kill $DIALOG_PID_OUT $DIALOG_PID_ERR 2>/dev/null
wait $DIALOG_PID_OUT $DIALOG_PID_ERR 2>/dev/null

# Mostrar resultados finales
dialog --title "Resultados Finales" --msgbox "Salida:\n$(tail -n 20 "$TEMP_OUT")\n\nErrores:\n$(cat "$TEMP_ERR")" 20 $OUT_WIDTH

# Limpiar
cleanup
