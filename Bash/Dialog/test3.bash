#!/bin/bash
var=$(dialog --stdout --menu "MENU" 10 30 3 1
red 2 green 3 blue)
case $? in
0) echo "Escogistes la opcion de aceptar"
case $var in
1) dialog --msgbox "Escogistes la
opcion 1: RED" 5 35;;
2) dialog --msgbox "Escogistes la
opcion 2: GREEN" 5 35;;
3) dialog --msgbox "Escogistes la
opcion 3: BLUE" 5 35;;
esac
;;
1)
dialog --sleep 4 --infobox "Escogistes
la opcion de cancelar" 8 45
dialog --pause "Saliendo del Sistema"
10 30 5
;;
esac
