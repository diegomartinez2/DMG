#!/usr/bin/env bash
# while-menu-dialog: a menu driven system information program

DIALOG_CANCEL=1
DIALOG_ESC=255
HEIGHT=0
WIDTH=0
# ------------------------------------------------------------------------------
display_load_file() {
FILENAME=$(dialog --stdout --backtitle "The Stochastic Self-Consistent Harmonic Approximation (SSCHA)" \
--title "Please choose a file" --fselect $HOME/ 14 48)
# echo "${FILENAME} file chosen."
}

display_result() {
  dialog --title "$1" \
    --no-collapse \
    --msgbox "$result" 0 0
}

display_help() {
  dialog --backtitle "The Stochastic Self-Consistent Harmonic Approximation (SSCHA)" \
  --title "SSCHA Help" --no-collapse --textbox sscha_help.txt 0 0
}

display_edit_file() {
  dialog --backtitle "The Stochastic Self-Consistent Harmonic Approximation (SSCHA)" \
  --title "Edit file content" --editbox sscha_help.txt 30 70 2> "${INPUT}"
}

# display_save_file() {
#   dialog --backtitle "The Stochastic Self-Consistent Harmonic Approximation (SSCHA)" \
#   --title "Save file" \
# --yesno {"Are you sure you want to overwrite " + $FILENAME + "?"} 7 60
#
# # Get exit status
# # 0 means user hit [yes] button.
# # 1 means user hit [no] button.
# # 255 means user hit [Esc] key.
# response=$?
# case $response in
#    0) cp ${INPUT} $FILENAME;;
#    1) echo "File not deleted.";;
#    255) echo "[ESC] key pressed.";;
# esac
# }

while true; do
  exec 3>&1
  selection=$(dialog \
    --backtitle "The Stochastic Self-Consistent Harmonic Approximation (SSCHA)" \
    --title "Menu" \
    --clear \
    --cancel-label "Exit" \
    --menu "Please select:" $HEIGHT $WIDTH 6 \
    "1" "Help" \
    "2" "Load SSCHA input file" \
    "3" "Edit input file" \
    "4" "Run calculation" \
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
      display_help
      ;;
    2 )
      display_load_file
      ;;
    3 )
      display_edit_file
      ;;
    4)
      display_save_file
      ;;
    5)
        result=$(echo "Scha:\n $FILENAME";uptime)
        display_result "System Information"
      ;;
  esac
done

# ---escribe en fichero----
# var="text to append";
# destdir=/some/directory/path/filename
#
# if [ -f "$destdir" ]
# then
#     echo "$var" > "$destdir"
# fi

# working with arrays
# for (( i=0; i<=${#myarray[@]}; i++ )); do
#      echo "${myarray[$i]}"
# done

# program output to window
# program | dialog --programbox 30 100
# --prgbox command height width   <--- mejor esto.
