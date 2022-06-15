#!/usr/bin/env bash
# while-menu-dialog: a menu driven system information program

DIALOG_CANCEL=1
DIALOG_ESC=255
HEIGHT=0
WIDTH=0

display_result() {
  dialog --title "$1" \
    --no-collapse \
    --msgbox "$result" 0 0
}

display_help() {
  dialog --title "SSCHA Help" --no-collapse --textbox sscha_help.txt 0 0

}

while true; do
  exec 3>&1
  selection=$(dialog \
    --backtitle "The Stochastic Self-Consistent Harmonic Approximation (SSCHA)" \
    --title "Menu" \
    --clear \
    --cancel-label "Exit" \
    --menu "Please select:" $HEIGHT $WIDTH 6 \
    "1" "SCHA input" \
    "2" "Calculator parameters" \
    "3" "Relax parameters" \
    "4" "Utilities" \
    "5" "Help" \
    "6" "Run calculation" \
    2>&1 1>&3)
  exit_status=$?
  exec 3>&-
  case $exit_status in
    $DIALOG_CANCEL)
      clear
      echo "Program terminated."
      exit
      ;;
    $DIALOG_ESC)
      clear
      echo "Program aborted." >&2
      exit 1
      ;;
  esac
  case $selection in
    1 )
      result=$(echo "Hostname: $HOSTNAME"; uptime)
      display_result "System Information"
      ;;
    2 )
      result=$(df -h)
      display_result "Disk Space"
      ;;
    3 )
      if [[ $(id -u) -eq 0 ]]; then
        result=$(du -sh /home/* 2> /dev/null)
        display_result "Home Space Utilization (All Users)"
      else
        result=$(du -sh $HOME 2> /dev/null)
        display_result "Home Space Utilization ($USER)"
      fi
      ;;
    4)
        result=$(echo "Hostname: $HOSTNAME"; uptime)
        display_result "System Information"
      ;;
    5)
        result=$(echo "Hostname: $HOSTNAME"; uptime)
        display_result "System Information"
      ;;
    6)
        result=$(echo "Hostname: $HOSTNAME"; uptime)
        display_result "System Information"
      ;;
  esac
done
