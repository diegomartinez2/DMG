#!/bin/bash
# Basically this does automatically:
#pkill -STOP lmp_mpi
#pkill -CONT lmp_mpi
# Nombre del proceso a buscar
PROCESO="lmp_mpi"

# Obtener el PID del primer proceso que coincida
PID=$(pgrep -f "$PROCESO" | head -n 1)

# Verificar si el proceso existe
if [ -z "$PID" ]; then
    echo "Error: No se encontró ningún proceso que coincida con '$PROCESO'."
    exit 1
fi

# Obtener el estado del proceso
# El flag 'stat' devuelve una letra (R, S, D, T, Z...)
ESTADO=$(ps -o stat= -p "$PID" | cut -c1)

case $ESTADO in
    T)
        echo "Proceso $PROCESO detectado como PAUSADO (T). Reanudando..."
        pkill -CONT -f "$PROCESO"
        echo "Señal SIGCONT enviada."
        ;;
    R|S|D)
        echo "Proceso $PROCESO detectado como ACTIVO ($ESTADO). Pausando..."
        pkill -STOP -f "$PROCESO"
        echo "Señal SIGSTOP enviada."
        ;;
    *)
        echo "El proceso está en un estado no gestionable ($ESTADO) o es un Zombie."
        ;;
esac
