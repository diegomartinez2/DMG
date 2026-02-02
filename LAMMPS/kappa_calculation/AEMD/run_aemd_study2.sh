#!/bin/bash

###############################################################################
# SCRIPT DE AUTOMATIZACIÓN AEMD - ESTUDIO SISTEMÁTICO SIN SOBRESCRITURA
###############################################################################

# 1. LISTAS DE PARÁMETROS
# -----------------------------------------------------------------------------
SIZES=("4 2 2" "8 2 2" "4 4 2" "16 4 4")
DELTA_TS=(5 10 20)
TIMESTEPS=("0.5" "1.0")
REPLICAS=3

INPUT_BASE="Thermo_AEMD_CuBHT_working.in"
LAMMPS_EXEC="lmp_serial" # Cambiar a "mpirun -np X lmp_mpi" si es necesario

# Archivos de soporte necesarios en la carpeta raíz
DATA_FILE="AA-supercell-optimized.lmp"
PARAMS_FILE="params.txt"

# 2. BUCLES DE EJECUCIÓN
# -----------------------------------------------------------------------------
for SIZE in "${SIZES[@]}"; do
    # Formatear el nombre del tamaño (ej. 4 2 2 -> 4x2x2)
    S_LABEL=$(echo $SIZE | tr ' ' 'x')

    for DT in "${DELTA_TS[@]}"; do
        #T_HOT=$(echo "300 + $DT" | bc -l)
        #T_COLD=$(echo "300 - $DT" | bc -l)
        T_HOT=$((300 + DT_VAL))
        T_COLD=$((300 - DT_VAL))

        for TS in "${TIMESTEPS[@]}"; do

            for ((R=1; R<=REPLICAS; R++)); do

                # CREACIÓN DE DIRECTORIO ÚNICO (Jerarquía: Tamaño/DT/TS/Replica)
                # Esto evita que los archivos .dat o .log se mezclen
                CASE_DIR="RESULTS/SIZE_${S_LABEL}/DT_${DT}/TS_${TS}/REP_${R}"
                mkdir -p "$CASE_DIR"

                # Nombre de los archivos de salida específicos para este caso
                OUTPUT_DAT="delta_T_S${S_LABEL}_DT${DT}_TS${TS}_R${R}.dat"
                LOG_FILE="log_S${S_LABEL}_DT${DT}_TS${TS}_R${R}.lammps"

                # Semillas aleatorias únicas
                S1=$((RANDOM + R * 7))
                S2=$((RANDOM + R * 13))
                S3=$((RANDOM + R * 17))

                echo "----------------------------------------------------"
                echo "EJECUTANDO: $CASE_DIR"
                echo "----------------------------------------------------"

                # Generar el archivo de entrada modificado
                # 1. Cambiamos replicate, dt y temperaturas
                # 2. Cambiamos las semillas de velocity
                # 3. IMPORTANTE: Cambiamos el nombre del archivo de salida en 'fix ave/time'
                sed -e "s/^replicate.*/replicate $SIZE/" \
                    -e "s/^variable dt equal.*/variable dt equal $TS/" \
                    -e "s/^variable T_hot_pulse equal.*/variable T_hot_pulse equal $T_HOT/" \
                    -e "s/^variable T_cold_pulse equal.*/variable T_cold_pulse equal $T_COLD/" \
                    -e "s/velocity all create \${T_eq} [0-9]*/velocity all create \${T_eq} $S1/" \
                    -e "s/velocity hot create \${T_hot_pulse} [0-9]*/velocity hot create \${T_hot_pulse} $S2/" \
                    -e "s/velocity cold create \${T_cold_pulse} [0-9]*/velocity cold create \${T_cold_pulse} $S3/" \
                    -e "s/fix output all ave\/time \(.*\) file .*/fix output all ave\/time \1 file $OUTPUT_DAT/" \
                    "$INPUT_BASE" > "$CASE_DIR/input.in"

                # Link simbólico a archivos pesados para no duplicar espacio en disco
                ln -sf "$(pwd)/$DATA_FILE" "$CASE_DIR/$DATA_FILE"
                ln -sf "$(pwd)/$PARAMS_FILE" "$CASE_DIR/$PARAMS_FILE"

                # Ejecución
                cd "$CASE_DIR"
                $LAMMPS_EXEC -in input.in -log "$LOG_FILE"
                cd - > /dev/null

            done
        done
    done
done

echo "Estudio completado. Los resultados están en la carpeta 'RESULTS/'"
