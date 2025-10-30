#!/bin/bash

# Este script extrae el último bloque de datos delimitado por marcadores.
# El archivo de entrada debe ser pasado como el primer argumento ($1).

# 1. Verificar si se proporcionó un archivo de entrada.
if [ -z "$1" ]; then
    echo "Error: Debe proporcionar el nombre del archivo de LAMMPS como argumento."
    echo "Uso: $0 nombre_del_archivo.log"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="out1"

# Marcadores para la extracción. Se usa una expresión regular flexible
# para manejar la variación en el espaciado entre las columnas.
START_MARKER="Step[[:space:]]+c_tl"
END_MARKER="Loop time of"

# La estrategia es:
# 1. Invertir el archivo con 'tac'. Esto hace que el "último" bloque sea el "primero".
# 2. Usar 'sed' para buscar el patrón del bloque (ahora invertido).
# 3. El patrón de inicio del bloque invertido es el END_MARKER original ("Loop time of").
# 4. El patrón de fin del bloque invertido es el START_MARKER original ("Step...").
# 5. La opción '-n' y el comando 'p' de sed se usan para imprimir el bloque.
# 6. 'q' (quit) se usa para detener el proceso después de encontrar el primer bloque (que es el último en el archivo original).
# 7. Invertir la salida de nuevo con 'tac' para restaurar el orden de tiempo.
# 8. Eliminar el marcador de fin ("Loop time of") que ahora está en la última línea del bloque.

# NOTA: La línea "Step..." es inclusiva, la línea "Loop time of" debe ser exclusiva.

tac "$INPUT_FILE" | \
sed -n "/${END_MARKER}/,/^.*${START_MARKER}/ {
    /^Loop time of/ b skip_end
    p
    /^.*${START_MARKER}/ q
    :skip_end
}" | \
tac > "$OUTPUT_FILE"

# Verificación de éxito.
if [ -s "$OUTPUT_FILE" ]; then
    echo "✅ El último bloque de datos ha sido extraído y guardado en '$OUTPUT_FILE'."
else
    echo "❌ No se pudo encontrar el bloque de datos o el archivo de entrada está vacío."
fi

# Ejecutar gnuplot y capturar toda la salida (stdout + stderr)
echo "Ejecutando gnuplot series.gp..."
gnuplot series.gp > gnuplot.out 2>&1

# Verificar si gnuplot falló
if [[ $? -ne 0 ]]; then
    echo "Error: gnuplot falló. Revisa 'gnuplot.out' para más detalles." >&2
    exit 1
fi

echo "Gnuplot ejecutado correctamente. Salida guardada en 'gnuplot.out'"

# Opcional: mostrar resumen
lines_in_data=$(wc -l < out1)
echo "Resumen: $lines_in_data líneas de datos extraídas."
cat gnuplot.out

exit 0
