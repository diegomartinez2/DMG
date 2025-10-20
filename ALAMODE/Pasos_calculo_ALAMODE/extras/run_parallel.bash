#!/bin/bash

# =================================================================
# CONFIGURACIÓN (Ajustar según sea necesario)
# =================================================================
LMP_EXEC="./lmp_mpi_chimes"
INPUT_SW="in.sw"
TOOL_PATH="/home/dmartinezgutierrez/ALAMODE/alamode-master/tools/extract.py"
LAMMPS_DAT="AA-supercell-optimized.lmp"
MAX_JOBS=4  # Número MÁXIMO de trabajos (simulaciones) concurrentes

# =================================================================
# FUNCIÓN DE LIMPIEZA
# =================================================================

# Esta función se ejecuta automáticamente al salir o en caso de error (trap EXIT).
cleanup() {
    echo "Realizando limpieza de archivos temporales..."
    # Limpia archivos temporales de configuración (lmp_config_*)
    find /tmp -user $USER -name "lmp_config_*.lammps" -delete 2>/dev/null
    # Limpia archivos temporales de salida (XFSET_run_*)
    rm -f XFSET_run_* 2>/dev/null
}
trap cleanup EXIT

CURRENT_JOBS=0

# =================================================================
# BUCLE 1: CÁLCULOS ARMÓNICOS (25 iteraciones)
# =================================================================

echo "Iniciando cálculos armónicos en paralelo (Max $MAX_JOBS trabajos)..."

for ((i=1;i<=25;i++))
do
    num=$(printf "%02d" $i)
    CONFIG_NAME="harm${num}.lammps"
    FINAL_OUTPUT_NAME="XFSET.harm${num}"

    # 1. Archivo Temporal Único para la CONFIGURACIÓN (Input de datos)
    TEMP_CONFIG=$(mktemp --tmpdir lmp_config_XXXXXX.lammps)
    cp "$CONFIG_NAME" "$TEMP_CONFIG"

    # 2. Nombre Único para el archivo de SALIDA (XFSET)
    # Usamos el PID ($$) del sub-shell como parte del nombre para asegurar unicidad
    TEMP_OUTPUT="XFSET_run_${num}_$$"

    # 3. EJECUTAR EN PARALELO (Sub-shell y ampersand &)
    (
        # Leemos 'in.sw' y aplicamos DOS sustituciones 'sed' encadenadas:
        # a) Reemplazar el archivo de datos ('tmp.lammps')
        # b) Reemplazar el nombre del archivo de salida ('XFSET')

        # NOTA: Las barras en $TEMP_CONFIG se escapan para 'sed'
        sed -e "s/tmp\.lammps/${TEMP_CONFIG//\//\\/}/g" \
            -e "s/XFSET/${TEMP_OUTPUT}/g" \
            "$INPUT_SW" | "$LMP_EXEC"

        # Copiar la salida temporal ($TEMP_OUTPUT) al nombre final
        cp "$TEMP_OUTPUT" "$FINAL_OUTPUT_NAME"

        # Limpiar los archivos temporales específicos de esta tarea
        rm -f "$TEMP_CONFIG" "$TEMP_OUTPUT"

    ) & # Envía el trabajo a segundo plano

    # 4. GESTIÓN DEL LÍMITE DE TRABAJOS
    CURRENT_JOBS=$((CURRENT_JOBS + 1))
    if [ "$CURRENT_JOBS" -ge "$MAX_JOBS" ]; then
        wait -n # Espera a que termine *uno* de los trabajos en background
        CURRENT_JOBS=$((CURRENT_JOBS - 1))
    fi
done

# Esperar a que terminen los trabajos restantes del bucle armónico
wait
echo "Todos los cálculos armónicos han finalizado."

# =================================================================
# BUCLE 2: CÁLCULOS CÚBICOS (7570 iteraciones)
# =================================================================

CURRENT_JOBS=0
echo "Iniciando cálculos cúbicos en paralelo (Max $MAX_JOBS trabajos)..."

for ((i=1;i<=7570;i++))
do
    num=$(printf "%04d" $i)
    CONFIG_NAME="cubic${num}.lammps"
    FINAL_OUTPUT_NAME="XFSET.cubic${num}"

    # 1. Archivo Temporal Único para la CONFIGURACIÓN (Input de datos)
    TEMP_CONFIG=$(mktemp --tmpdir lmp_config_XXXXXX.lammps)
    cp "$CONFIG_NAME" "$TEMP_CONFIG"

    # 2. Nombre Único para el archivo de SALIDA (XFSET)
    TEMP_OUTPUT="XFSET_run_${num}_$$"

    # 3. EJECUTAR EN PARALELO
    (
        # Doble sustitución 'sed' para entrada y salida
        sed -e "s/tmp\.lammps/${TEMP_CONFIG//\//\\/}/g" \
            -e "s/XFSET/${TEMP_OUTPUT}/g" \
            "$INPUT_SW" | "$LMP_EXEC"

        # Copiar la salida temporal al nombre final
        cp "$TEMP_OUTPUT" "$FINAL_OUTPUT_NAME"

        # Limpiar los archivos temporales
        rm -f "$TEMP_CONFIG" "$TEMP_OUTPUT"

    ) & # Envía el trabajo a segundo plano

    # 4. GESTIÓN DEL LÍMITE DE TRABAJOS
    CURRENT_JOBS=$((CURRENT_JOBS + 1))
    if [ "$CURRENT_JOBS" -ge "$MAX_JOBS" ]; then
        wait -n # Espera a que termine *uno* de los trabajos en background
        CURRENT_JOBS=$((CURRENT_JOBS - 1))
    fi
done

# Esperar a que terminen los trabajos restantes del bucle cúbico
wait
echo "Todos los cálculos cúbicos han finalizado."

# =================================================================
# POST-PROCESAMIENTO
# =================================================================

echo "Iniciando post-procesamiento con Python..."

# Procesamiento armónico
python "$TOOL_PATH" --LAMMPS="$LAMMPS_DAT" XFSET.harm* > DFSET_harmonic

# Procesamiento cúbico
python "$TOOL_PATH" --LAMMPS="$LAMMPS_DAT" XFSET.cubic* > DFSET_cubic

# Fusión
cat DFSET_harmonic DFSET_cubic > DFSET_merged

echo "Proceso completado. Archivo final: DFSET_merged"
