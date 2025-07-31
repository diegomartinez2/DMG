# run_all_lammps.sh
#!/bin/bash

# Asegúrate de que LAMMPS esté en tu PATH o especifica la ruta completa
# Por ejemplo, si lo instalaste con Conda:
# LAMMPS_EXE=$(which lmp)
LAMMPS_EXE="lmp" # Asume que 'lmp' está en tu PATH

DATA_PREFIX="displaced_config_"
LAMMPS_IN="run_lammps_single.in"
OUTPUT_DIR="lammps_outputs" # Directorio para guardar los resultados de LAMMPS

mkdir -p $OUTPUT_DIR # Crea el directorio si no existe

# Obtén el número de archivos de configuración generados por el script de Python
NUM_CONFIGS=$(ls ${DATA_PREFIX}*.lammps | wc -l)

if [ "$NUM_CONFIGS" -eq 0 ]; then
    echo "Error: No se encontraron archivos de configuración desplazados (${DATA_PREFIX}*.lammps)."
    echo "Asegúrate de haber ejecutado 'generate_lammps_inputs.py' correctamente."
    exit 1
fi

echo "Procesando $NUM_CONFIGS configuraciones con LAMMPS..."

for i in $(seq 1 $NUM_CONFIGS); do
    # Formatear el número con ceros iniciales (ej., 0001, 0002)
    config_num=$(printf "%04d" $i)

    input_file="${DATA_PREFIX}${config_num}.lammps"
    output_force_file="${OUTPUT_DIR}/forces_${config_num}.dat"

    echo "Calculando fuerzas para ${input_file}..."

    # Ejecuta LAMMPS, pasando las variables para el archivo de entrada y salida
    $LAMMPS_EXE -in ${LAMMPS_IN} -var structure_file ${input_file} -var output_file ${output_force_file} > ${OUTPUT_DIR}/lammps_log_${config_num}.log

    if [ $? -eq 0 ]; then
        echo "Fuerzas guardadas en ${output_force_file}"
    else
        echo "Error en el cálculo de LAMMPS para ${input_file}. Revisa el log."
    fi
done

echo "Todos los cálculos de LAMMPS completados."
