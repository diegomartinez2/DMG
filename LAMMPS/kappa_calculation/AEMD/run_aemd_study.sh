#!/bin/bash

###############################################################################
# SCRIPT DE AUTOMATIZACIÓN PARA ESTUDIO SISTEMÁTICO CuBHT AEMD (LAMMPS)
###############################################################################

# 1. LISTAS DE PARÁMETROS (Edita estas listas según tus necesidades)
# -----------------------------------------------------------------------------
# Tamaños de super-red: "X Y Z"
SIZES=("4 2 2" "8 2 2" "4 4 2" "4 2 4" "4 4 4" "8 4 4" "16 2 2" "16 4 4")

# Saltos de temperatura (Delta T alrededor de T_eq=300K)
# Ejemplo: 10 significa T_hot=310 y T_cold=290
DELTA_TS=(5 15 25)

# Pasos de tiempo (Time-steps en fs)
TIMESTEPS=(0.5 1.0 1.5)

# Número de réplicas por cada configuración (cambia las semillas aleatorias)
REPLICAS=3

# Archivo de entrada base y ejecutables
INPUT_BASE="Thermo_AEMD_CuBHT_working.in"
LAMMPS_EXEC="lmp_serial" # O mpirun -np 4 lmp_mpi, etc.

# 2. BUCLES DE EJECUCIÓN
# -----------------------------------------------------------------------------

for SIZE in "${SIZES[@]}"; do
    # Convertir "4 2 2" a "422" para nombres de carpeta
    SIZE_NAME=$(echo $SIZE | tr -d ' ')

    for DT_VAL in "${DELTA_TS[@]}"; do
        #T_HOT=$(echo "300 + $DT_VAL" | bc -l)
        #T_COLD=$(echo "300 - $DT_VAL" | bc -l)
        #T_HOT=$(bc -l <<< "300 + $DT_VAL")
        #T_COLD=$(bc -l <<< "300 - $DT_VAL")
        # Opción limpia con aritmética bash (recomendada aquí)
        T_HOT=$((300 + DT_VAL))
        T_COLD=$((300 - DT_VAL))

        # ────────────────────────────────────────────────
        # Imprimir para depurar / verificar
        echo "  ΔT = $DT_VAL  →  T_hot = $T_HOT K    T_cold = $T_COLD K"
        # ────────────────────────────────────────────────

        for STEP in "${TIMESTEPS[@]}"; do

            for ((R=1; R<=REPLICAS; R++)); do

                # Crear identificador único y carpeta de trabajo
                WORK_DIR="run_S${SIZE_NAME}_DT${DT_VAL}_TS${STEP}_R${R}"
                echo "==> Iniciando: $WORK_DIR"

                mkdir -p "$WORK_DIR"

                # Generar semillas aleatorias únicas para esta réplica
                SEED1=$((RANDOM + R * 100))
                SEED2=$((RANDOM + R * 200))
                SEED3=$((RANDOM + R * 300))

                # Crear el archivo de entrada modificado para esta ejecución
                # Usamos sed para reemplazar las líneas específicas del script original
                sed -e "s/^replicate.*/replicate $SIZE/" \
                    -e "s/^variable dt equal.*/variable dt equal $STEP/" \
                    -e "s/^variable T_hot_pulse equal.*/variable T_hot_pulse equal $T_HOT/" \
                    -e "s/^variable T_cold_pulse equal.*/variable T_cold_pulse equal $T_COLD/" \
                    -e "s/velocity all create \${T_eq} [0-9]*/velocity all create \${T_eq} $SEED1/" \
                    -e "s/velocity hot create \${T_hot_pulse} [0-9]*/velocity hot create \${T_hot_pulse} $SEED2/" \
                    -e "s/velocity cold create \${T_cold_pulse} [0-9]*/velocity cold create \${T_cold_pulse} $SEED3/" \
                    "$INPUT_BASE" > "$WORK_DIR/input.in"

                # Copiar archivos necesarios (potenciales, datos, etc.)
                cp params.txt AA-supercell-optimized.lmp "$WORK_DIR/"

                # Ejecutar LAMMPS
                cd "$WORK_DIR"
                $LAMMPS_EXEC -in input.in > lammps.log

                # Opcional: Limpiar archivos pesados si la simulación termina bien
                # rm params.txt AA-supercell-optimized.lmp

                cd ..

                echo "==> Finalizado: $WORK_DIR"
            done
        done
    done
done

echo "Estudio sistemático completado."
