#!/bin/bash

# --- Función de Extracción (Proporcionada por el usuario) ---
extract() {
    if [ -f "$1" ] ; then
        case "$1" in
            *.tar.bz2)   tar xvjf "$1"     ;;
            *.tar.gz)    tar xvzf "$1"     ;;
            *.bz2)       bunzip2 "$1"      ;;
            *.rar)       unrar x "$1"      ;;
            *.gz)        gunzip "$1"       ;;
            *.tar)       tar xvf "$1"      ;;
            *.tbz2)      tar xvjf "$1"     ;;
            *.tgz)       tar xvzf "$1"     ;;
            *.zip)       unzip "$1"        ;;
            *.Z)         uncompress "$1"   ;;
            *.7z)        7z x "$1"         ;;
            *)           echo "'$1' no se puede extraer con >extract<" ;;
        esac
    else
        echo "'$1' ¡no es un archivo válido!"
    fi
}

# --- Función de Compresión ---
# Recibe: $1 (archivo origen), $2 (archivo destino/salida)
compress() {
    local source=$1
    local destination=$2

    if [ ! -e "$source" ]; then
        echo "Error: El archivo o carpeta '$source' no existe."
        exit 1
    fi

    echo "Comprimiendo '$source' en '$destination'..."

    case "$destination" in
        *.tar.bz2|*.tbz2) tar cvjf "$destination" "$source" ;;
        *.tar.gz|*.tgz)   tar cvzf "$destination" "$source" ;;
        *.tar)            tar cvf "$destination" "$source"  ;;
        *.zip)            zip -r "$destination" "$source"   ;;
        *.7z)             7z a "$destination" "$source"     ;;
        *.gz)             gzip -c "$source" > "$destination" ;;
        *.bz2)            bzip2 -c "$source" > "$destination" ;;
        *)                echo "Formato de compresión '$destination' no soportado." ;;
    esac
}

# --- Lógica Principal del Script ---

if [ $# -eq 1 ]; then
    # Un solo argumento: Descomprimir
    extract "$1"
elif [ $# -eq 2 ]; then
    # Dos argumentos: Comprimir (Origen Destino)
    compress "$1" "$2"
else
    # Uso incorrecto
    echo "Uso del script:"
    echo "  Descomprimir: $0 archivo.ext"
    echo "  Comprimir:    $0 origen destino.ext"
    echo ""
    echo "Ejemplos:"
    echo "  $0 mis_datos.zip             (Extrae el contenido)"
    echo "  $0 carpeta/ backup.tar.gz    (Comprime la carpeta)"
    exit 1
fi
