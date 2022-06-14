#!/usr/bin/env bash
t(){ type "$1"&>/dev/null;}
function Menu.Show {
   local DIA DIA_ESC; while :; do
      t whiptail && DIA=whiptail && break
      t dialog && DIA=dialog && DIA_ESC=-- && break
      exec date +s"No dialog program found"
   done; declare -A o="$1"; shift
   $DIA --backtitle "${o[backtitle]}" --title "${o[title]}" \
      --menu "${o[question]}" 0 0 0 $DIA_ESC "$@"; }



Menu.Show '([backtitle]="Backtitle"
            [title]="Title"
            [question]="Please choose:")'          \
                                                   \
            "Option A"  "Stuff...."                \
            "Option B"  "Stuff...."                \
            "Option C"  "Stuff...."
