#!/usr/bin/env bash
# creamos la funcion fundialog
fundialog=${fundialog=dialog}
fecha=`$fundialog --stdout -- backtitle "SSCHA press F1 for help"--title "Welcome" --msgbox "Welcome to the SSCHA code." 0 0`
echo $fecha
