#!/bin/bash
# Por Ignacio Alba Obaya
# https://aplicacionesysistemas.com

# Creamos la varaible funcheck en la que almacenamos la
# orden dialog con la opción --separate-output
funcheck=(dialog --separate-output --checklist "Selecciona los grupos a los que pertenece:" 0 0 0)

# Definimos las opciones que apareceran en pantalla
# aparecerán encendidas las que marquemos con on.
opciones=(1 "opción 1" on
 2 "opción 2" off
 3 "opción 3" off
 4 "opción 4" off
 5 "opción 5" on
 6 "opción 6" off
 7 "opción 7" off)

# Creamos la funcion selecciones que ejecuta funcheck con opciones
# y reenvia la salida al terminal para que el for siguiente ejecute
# los comandos
selecciones=$("${funcheck[@]}" "${opciones[@]}" 2>&1 >/dev/tty)

# limpiamos la pantalla
clear

# añadimos un for para que ejecute un comando en función de
# las selecciones realizadas puedes cambiar el echo por
# cualquier comando o secuencias de comandos
for seleccion in $selecciones
do
 case $seleccion in
 1)
 echo "Escogiste la opción 1"
 ;;
 2)
 echo "Escogiste la opción 2"
 ;;
 3)
 echo "Escogiste la opción 3"
 ;;
 4)
 echo "Escogiste la opción 4"
 ;;
 5)
 echo "Escogiste la opción 5"
 ;;
 6)
 echo "Escogiste la opción 6"
 ;;
 7)
 echo "Escogiste la opción 7"
 ;;
 esac
done
