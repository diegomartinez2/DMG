#! /bin/sh
a="Mr."
b="Tux"
c="Penguin"
dialog --title "Title" --backtitle "Customer Data" --ok-label "Save" \
  --stdout --form "Catalog" 10 60 3 \
  "Salutation  " 1 1 "$a" 1 15 30 0 \
  "Family Name " 2 1 "$b" 2 15 30 0 \
  "First Name  " 3 1 "$c" 3 15 30 0 > output.txt
a=$(cat output.txt | head -1)
b=$(cat output.txt | head -2 | tail -1)
c=$(cat output.txt | head -3 | tail -1)
rm output.txt
dialog --title "Title" --backtitle "Background Title" --msgbox \
  "Saved values: \n $a \n $b \n $c " 0 0
