#!/bin/bash

# Archivos temporales para salida y errores
TEMP_OUT=$(mktemp)
TEMP_ERR=$(mktemp)

# Función para limpiar archivos temporales
cleanup() {
    rm -f "$TEMP_OUT" "$TEMP_ERR"
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

# Función para monitorear salida y errores
monitor_output() {
    # Ventana para salida estándar (90% de la altura)
    dialog --title "Salida Estándar - LAMMPS" --tailboxbg "$TEMP_OUT" $OUT_HEIGHT $OUT_WIDTH \
        --begin 2 2 &

    DIALOG_PID_OUT=$!

    # Ventana para errores (debajo, más pequeña)
    dialog --title "Errores" --tailboxbg "$TEMP_ERR" $ERR_HEIGHT $ERR_WIDTH \
        --begin $((OUT_HEIGHT + SEP_HEIGHT + 2)) 2 &

    DIALOG_PID_ERR=$!
}

# Función para ejecutar el cálculo LAMMPS con barra de progreso
run_lammps() {
    # Suponiendo que tienes un archivo de entrada 'in.lammps'
    # y que LAMMPS genera un archivo de log (ej. log.lammps)
    LAMMPS_INPUT="in.lammps"
    LAMMPS_LOG="log.lammps"

    # Obtener número total de pasos del archivo de entrada (ajusta según tu input)
    TOTAL_STEPS=$(grep -i "run" "$LAMMPS_INPUT" | grep -o '[0-9]\+' | head -1)
    [ -z "$TOTAL_STEPS" ] && TOTAL_STEPS=100  # Valor por defecto si no se encuentra

    # Ejecutar LAMMPS y monitorear progreso
    {
        lammps_command="lmp_serial -in $LAMMPS_INPUT"  # Cambia por tu comando LAMMPS (ej. mpirun -np 4 lmp_mpi)
        $lammps_command > "$TEMP_OUT" 2> >(while read -r line; do echo "$line" >> "$TEMP_ERR"; done) &

        LAMMPS_PID=$!

        # Monitorear el archivo de log para el progreso
        while [ -d /proc/$LAMMPS_PID ]; do
            if [ -f "$LAMMPS_LOG" ]; then
                # Extraer el paso actual (ajusta según el formato de tu log)
                CURRENT_STEP=$(tail -n 10 "$LAMMPS_LOG" | grep -oE '^[0-9]+' | tail -n 1)
                [ -z "$CURRENT_STEP" ] && CURRENT_STEP=0
                # Calcular porcentaje
                PERCENT=$((CURRENT_STEP * 100 / TOTAL_STEPS))
                [ $PERCENT -gt 100 ] && PERCENT=100
                echo $PERCENT
            else
                echo 0
            fi
            sleep 1
        done
    } | dialog --title "Progreso del Cálculo LAMMPS" --gauge "Ejecutando simulación..." 8 $OUT_WIDTH 0

    # Esperar a que LAMMPS termine
    wait $LAMMPS_PID

    # Mostrar resultados finales
    dialog --title "Resultados Finales" --msgbox "Salida:\n$(tail -n 20 "$TEMP_OUT")\n\nErrores:\n$(cat "$TEMP_ERR")" 20 $OUT_WIDTH
}

# Iniciar monitoreo de salida y errores
monitor_output

# Ejecutar LAMMPS
run_lammps

# Terminar procesos de dialog
kill $DIALOG_PID_OUT $DIALOG_PID_ERR 2>/dev/null
wait $DIALOG_PID_OUT $DIALOG_PID_ERR 2>/dev/null

# Limpiar
cleanup
