#!/bin/bash
# Por Ignacio Alba Obaya
# https://aplicacionesysistemas.com
# Almacenar el resultado de dialog en una variable
# creamos la funcion fundialog
fundialog=${fundialog=dialog}
# creamos una variable fecha con la salida de dialog
# redireccinando con --stdout la salida de dialog
# a la salida estandar, fijate en que la función esta entre
# acentos de los de la tecla [ si no es asi no funciona.
fecha=`$fundialog --stdout --title "Calendario" --calendar "Escoge una fecha" 0 0`
# Mostramos la fecha capturada
echo $fecha
