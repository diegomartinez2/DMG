#!/bin/bash
onoff=on
cmd=(dialog --output-fd 1 --separate-output --extra-button --extra-label 'Select All' --cancel-label 'Select None' --checklist 'Choose the tools to install:' 0 0 0)
load-dialog () {
    options=(
                1 'Option 1' $onoff
                2 'Option 2' $onoff
                3 'Option 3' $onoff
                4 'Option 4' $onoff
                5 'Option 5' $onoff
                6 'Option 6 (Depends Option 5)' $onoff
                7 'Option 7' $onoff
                8 'Option 8' $onoff
    )
    choices=$("${cmd[@]}" "${options[@]}")
}
load-dialog
exit_code="$?"
while [[ $exit_code -ne 0 ]]; do
case $exit_code in
    1) clear; onoff=off; load-dialog;;
    3) clear; onoff=on; load-dialog;;
esac
exit_code="$?"
done
clear
for choice in $choices
do
    case $choice in
        1) echo 'First Option';;
        2) echo 'Second Option';;
        3) echo 'Third Option';;
        4) echo 'Fourth Option';;
        5) echo 'Fifth Option';;
        6) echo 'Sixth Option';;
        7) echo 'Seventh Option';;
        8) echo 'Eighth Option';;
    esac
done
sleep infinity
