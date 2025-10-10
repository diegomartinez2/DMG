#!/bin/bash

# =======================================================
# 1. Configuración del Entorno
# =======================================================

# Ruta al ejecutable de LAMMPS (ajustar según tu sistema)
# Ejemplo para ejecución serial (lmp) o paralela (mpirun)
# LMP_EXE="lmp"
LMP_EXE="mpirun -np 4 lmp_mpi"

# Nombre del archivo de entrada de LAMMPS
LAMMPS_INPUT="kappa.in"

# Nombre del script de post-procesamiento de Python
PYTHON_SCRIPT="post_process.py"


# =======================================================
# 2. Limpieza y Ejecución de LAMMPS
# =======================================================

echo "Iniciando el cálculo de Conductividad Térmica..."
echo "Limpiando archivos de salida anteriores..."

# Eliminar archivos de salida antiguos
rm -f J0Jt.dat params.txt log.lammps

# Ejecutar LAMMPS
# El '-in' especifica el archivo de entrada.
echo "Ejecutando LAMMPS con $LMP_EXE..."
$LMP_EXE -in $LAMMPS_INPUT

# Verificar si LAMMPS se ejecutó correctamente
if [ $? -ne 0 ]; then
    echo "❌ Error: La simulación de LAMMPS falló. Revisar log.lammps."
    exit 1
fi

echo "✅ Simulación de LAMMPS completada."

# =======================================================
# 3. Post-Procesamiento con Python
# =======================================================

echo "Iniciando el post-procesamiento con Python..."

# Ejecutar el script de Python
python3 $PYTHON_SCRIPT

# Verificar si Python se ejecutó correctamente
if [ $? -ne 0 ]; then
    echo "❌ Error: El script de Python falló. Revisar post_process.py."
    exit 1
fi

echo "✅ Post-procesamiento completado. Resultados mostrados arriba."
