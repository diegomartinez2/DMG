#!/bin/bash

echo "Iniciando ejecución PARALELA..."

# El loop va de 1 a 10
for i in {1..10}; do
    # Formateo del número con cero inicial (01, 02, ..., 10)
    script_num=$(printf "%02d" $i)
    script_name="a${script_num}.py"

    echo "Lanzando al background: $script_name"

    # Ejecuta el programa y lo envía al background con &
    python "$script_name" &
done

# Esperar a que todos los procesos lanzados por el bucle terminen
echo "Esperando a que todos los 10 procesos finalicen..."
wait

echo "Ejecución paralela completada."
