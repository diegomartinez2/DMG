#!/bin/bash
# Por Ignacio Alba Obaya
# https://aplicacionesysistemas.com
# ejecutamos el cuadro de dialogo poniendo al final 2>/tmp/nombre.tmp.$$
# almacenamos en un archivo el nombre introducido.
# recordar que 2> redirecciona la salida de error a un archivo.
dialog --title "Nombre" --inputbox "Pon tu nombre:" 0 0 2>/tmp/nombre.tmp.$$
# borramos la pantalla
clear
# mostramos el nombre almacenado
cat /tmp/nombre.tmp.$$
# borramos el archivo con el nombre
rm -f /tmp/nombre.tmp.$$
# hacemos un salto de linea para que no se nos amontone con el prompt
echo -e "\n"
