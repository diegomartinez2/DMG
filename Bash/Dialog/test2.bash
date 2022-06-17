#!/bin/bash
var=$( dialog --stdout --checklist "Elija su color
preferido" 10 40 3 1 red on 2 green off 3 blue
on 4 pink off )
dialog --msgbox "Valores Capturados: $var" 5 50

# -------otra forma-----------
var=$( dialog --stdout --separate-output --
checklist "Elija su color preferido" 10 40 3 1 red
on 2 green off 3 blue on 4 pink off )
echo "Valores Capturados: "$var
