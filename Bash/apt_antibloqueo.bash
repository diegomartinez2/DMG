#!/bin/bash

echo "🔍 Verificando si packagekitd está activo..."
PK_PID=$(pgrep packagekitd)

if [ -n "$PK_PID" ]; then
    echo "⚠️ packagekitd está corriendo con PID $PK_PID. Deteniendo temporalmente..."
    sudo systemctl stop packagekit
    echo "✅ packagekit detenido."
else
    echo "✅ packagekitd no está corriendo."
fi

echo "🔄 Ejecutando 'sudo apt update'..."
sudo apt update

if [ -n "$PK_PID" ]; then
    echo "🔁 Reiniciando packagekit..."
    sudo systemctl start packagekit
    echo "✅ packagekit reiniciado."
fi

echo "🎉 Proceso completado."
