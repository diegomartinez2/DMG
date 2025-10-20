#!/bin/bash

# =================================================================
# Variables de Entorno (Ajusta si es necesario)
# =================================================================
LMP_EXEC="./lmp_mpi_chimes"
INPUT_SW="in.sw"
TOOL_PATH="/home/dmartinezgutierrez/ALAMODE/alamode-master/tools/extract.py"
LAMMPS_DAT="AA-supercell-optimized.lmp"


# Función de limpieza (importante en caso de errores)
cleanup() {
    # Eliminamos el archivo temporal principal y el input temporal de LAMMPS (si quedan)
    find /tmp -user $USER -name "lammps_input_*.tmp" -delete 2>/dev/null
    find /tmp -user $USER -name "lmp_config_*.lammps" -delete 2>/dev/null
}
# Trap para asegurar que la limpieza se ejecute al salir (incluso con error)
trap cleanup EXIT


# =================================================================
# PRIMER BUCLE: Archivos armónicos
# =================================================================

echo "Iniciando cálculos armónicos..."
# # Establecer límite de trabajos concurrentes (ejemplo: 4 CPUs)
# MAX_JOBS=4
# CURRENT_JOBS=0
#
# echo "Iniciando cálculos armónicos en paralelo (Max $MAX_JOBS trabajos)..."


for ((i=1;i<=25;i++))
do
    num=$(printf "%02d" $i)
    CONFIG_NAME="harm${num}.lammps"
    OUTPUT_NAME="XFSET.harm${num}"

    # 1. Crear un nombre de archivo temporal único para la configuración de LAMMPS
    # Usaremos una extensión real para la seguridad del entorno
    TEMP_CONFIG=$(mktemp --tmpdir lmp_config_XXXXXX.lammps)

    # 2. Copiar el archivo de configuración al nuevo nombre temporal
    cp "$CONFIG_NAME" "$TEMP_CONFIG"

    # 3. Leer el input 'in.sw', modificar el nombre del archivo de lectura,
    # y pasar el input modificado directamente a lammps.

    # Usamos sed para reemplazar 'tmp.lammps' (o el nombre original) por la ruta única.
    # El patrón de reemplazo debe ser robusto, asumiendo que el archivo de lectura
    # está en el in.sw, por ejemplo: "read_data tmp.lammps"

    echo "  -> Ejecutando $CONFIG_NAME con archivo temporal $TEMP_CONFIG"

    sed "s/tmp\.lammps/${TEMP_CONFIG//\//\\/}/g" "$INPUT_SW" | "$LMP_EXEC"

    # NOTA: La parte 's/tmp\.lammps/${TEMP_CONFIG//\//\\/}/g' se asegura
    # de que si $TEMP_CONFIG tiene barras diagonales (/), estas se escapen (\/) para sed.

    # 4. Copiar el archivo de salida
    cp XFSET "$OUTPUT_NAME"

    # 5. Borrar el archivo temporal de configuración inmediatamente después de usarlo
    rm -f "$TEMP_CONFIG"
done

# ... (El segundo bucle para 'cubic' seguiría una lógica idéntica) ...
# =================================================================
# SEGUNDO BUCLE: Archivos cúbicos (Esquema idéntico)
# =================================================================
echo "Iniciando cálculos cúbicos..."

for ((i=1;i<=7570;i++))
do
    num=$(printf "%04d" $i)
    CONFIG_NAME="cubic${num}.lammps"
    OUTPUT_NAME="XFSET.cubic${num}"

    TEMP_CONFIG=$(mktemp --tmpdir lmp_config_XXXXXX.lammps)
    cp "$CONFIG_NAME" "$TEMP_CONFIG"

    echo "  -> Ejecutando $CONFIG_NAME con archivo temporal $TEMP_CONFIG"

    sed "s/tmp\.lammps/${TEMP_CONFIG//\//\\/}/g" "$INPUT_SW" | "$LMP_EXEC"

    cp XFSET "$OUTPUT_NAME"
    rm -f "$TEMP_CONFIG"
done


# =================================================================
# POST-PROCESAMIENTO (Esta parte no necesita cambios)
# =================================================================
echo "Iniciando post-procesamiento con Python..."

# Procesamiento armónico
echo "python extract.py --LAMMPS=$LAMMPS_DAT XFSET.harm* > DFSET_harmonic"
python "$TOOL_PATH" --LAMMPS="$LAMMPS_DAT" XFSET.harm* > DFSET_harmonic

# Procesamiento cúbico
echo "python extract.py --LAMMPS=$LAMMPS_DAT XFSET.cubic* > DFSET_cubic"
python "$TOOL_PATH" --LAMMPS="$LAMMPS_DAT" XFSET.cubic* > DFSET_cubic

# Fusión
cat DFSET_harmonic DFSET_cubic > DFSET_merged

echo "Proceso completado. Archivo final: DFSET_merged"
