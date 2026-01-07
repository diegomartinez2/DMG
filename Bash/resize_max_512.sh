#!/bin/bash

# Script para redimensionar imágenes a un máximo de 512px en el lado más largo
# Manteniendo la relación de aspecto
# Requiere ImageMagick instalado

# Uso: ./resize_max_512.sh [directorio] [opcional: --backup]
# Si no se pasa directorio, usa el actual

DIRECTORIO="${1:-.}"
BACKUP=false

# Opción para hacer backup de originales
if [[ "$2" == "--backup" || "$1" == "--backup" ]]; then
    BACKUP=true
    [[ "$1" == "--backup" ]] && DIRECTORIO="${2:-.}"
fi

# Extensiones de imagen comunes (puedes añadir más)
EXTENSIONES="jpg jpeg png gif webp bmp tiff tif"

echo "Redimensionando imágenes en: $DIRECTORIO"
echo "Tamaño máximo: 512px (lado más largo)"
echo "Manteniendo aspect ratio"
[[ "$BACKUP" == true ]] && echo "Se creará backup de originales"

# Contador
contador=0

# Procesar cada extensión
for ext in $EXTENSIONES; do
    # Buscar archivos (case insensitive) y procesarlos
    find "$DIRECTORIO" -type f -iname "*.$ext" | while read -r imagen; do
        # Obtener dimensiones actuales
        ancho=$(magick identify -format "%w" "$imagen" 2>/dev/null)
        alto=$(magick identify -format "%h" "$imagen" 2>/dev/null)

        # Si no se pudo leer, saltar
        [[ -z "$ancho" || -z "$alto" ]] && continue

        lado_mas_largo=$(( ancho > alto ? ancho : alto ))

        # Solo redimensionar si el lado más largo supera 512px
        if [[ $lado_mas_largo -gt 512 ]]; then
            echo "Redimensionando: $imagen (${ancho}x${alto} → máx 512)"

            if [[ "$BACKUP" == true ]]; then
                cp "$imagen" "$imagen.bak"
            fi

            # Redimensionar: -resize 512x512> significa "máximo 512 en cualquier eje"
            magick "$imagen" -resize 512x512\> "$imagen"

            ((contador++))
        fi
    done
done

echo "¡Listo! Se redimensionaron $contador imágenes."
[[ "$BACKUP" == true ]] && echo "Backups guardados con extensión .bak"
