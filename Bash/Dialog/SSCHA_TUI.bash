#!/usr/bin/env bash
# while-menu-dialog: a menu driven system information program

DIALOG_CANCEL=1
DIALOG_ESC=255
HEIGHT=0
WIDTH=0
# ------------------------------------------------------------------------------
# scha input

__SCHA_LAMBDA_A__=0.5
__SCHA_LAMBDA_W__=0.5
__SCHA_MINSTRUC__=.false.
__SCHA_PRECOND_WYCK__=.true.
__SCHA_PRECOND_DYN__=.true.
__SCHA_ROOTREP__="normal"
__SCHA_NEGLECT_SYMMETRIES__=.false.
__SCHA_NRANDOM_EFF__=500
__SCHA_NRANDOM__=1000
__SCHA_MEANINGFUL__=1e-4
__SCHA_EQENERGY__=-144.40680397
__SCHA_FILDYN__="../ensemble_data_test/dyn"
__SCHA_NQIRR__=1
__SCHA_DATADIR__="../ensemble_data_test"
__SCHA_ISBIN__=
__SCHA_T__=0.0d0
__SCHA_TG__=0
__SCHA_SUPERCELLSIZE__=" 1 1 1 "
__SCHA_MAXSTEPS__=80
__SCHA_STRESSOFFSET__=
__SCHA_GRADIOP__="all"
__SCHA_POPULATION__=1
__SCHA_PRINTSTRESS__=.true.
__SCHA_USESPGLIB__=.false.
display_Scha_input() {
  # open fd
  exec 3>&1

  # Store data to $SCHA_VALUES variable  later put into array with myarray=($myvar)  or  read -a myarray <<< $myvar
  SCHA_VALUES=$(dialog --ok-label "Submit" \
	  --backtitle "The Stochastic Self-Consistent Harmonic Approximation (SSCHA)" \
	  --title "SCHA parameters" \
	  --form "Fildyn, nqirr and T are mandatory." \
    30 50 0 \
     "fildyn_prefix:"      1 1	"$__SCHA_FILDYN__"	                 1 20 30 0 \
     "nqirr:"              2 1  "$__SCHA_NQIRR__"                    2 20 30 0 \
     "T:"                  3 1  "$__SCHA_T__"                        3 20 30 0 \
     "lambda_a:"           4 1  "$__SCHA_LAMBDA_A__"                 4 20 30 0 \
     "lambda_w:"           5 1  "$__SCHA_LAMBDA_W__"                 5 20 30 0 \
     "minim_struc:"        6 1  "$__SCHA_MINSTRUC__"                 6 20 30 0 \
     "precond_wyck:"       7 1  "$__SCHA_PRECOND_WYCK__"             7 20 30 0 \
     "preconditioning:"    8 1  "$__SCHA_PRECOND_DYN__"              8 20 30 0 \
     "root_representation:" 9 1  "$__SCHA_ROOTREP__"                  9 20 30 0 \
     "neglect_symmetries:" 10 1  "$__SCHA_NEGLECT_SYMMETRIES__"      10 20 30 0 \
     "n_random_eff:"       11 1  "$__SCHA_NRANDOM_EFF__"             11 20 30 0 \
     "n_random:"           12 1  "$__SCHA_NRANDOM__"                 12 20 30 0 \
     "meaningful_factor:"  13 1  "$__SCHA_MEANINGFUL__"              13 20 30 0 \
     "eq_energy:"          14 1  "$__SCHA_EQENERGY__"                14 20 30 0 \
     "data_dir:"           15 1  "$__SCHA_DATADIR__"                 15 20 30 0 \
     "load_bin:"           16 1  "$__SCHA_ISBIN__"                   16 20 30 0 \
     "Tg:"                 17 1  "$__SCHA_TG__"                      17 20 30 0 \
     "supercell_size:"     18 1  "$__SCHA_SUPERCELLSIZE__"           18 20 30 0 \
     "max_ka:"             19 1  "$__SCHA_MAXSTEPS__"                19 20 30 0 \
     "stress_offset:"      20 1  "$__SCHA_STRESSOFFSET__"            20 20 30 0 \
     "gradi_op:"           21 1  "$__SCHA_GRADIOP__"                 21 20 30 0 \
     "population:"         22 1  "$__SCHA_POPULATION__"              22 20 30 0 \
     "print_stress:"       23 1  "$__SCHA_PRINTSTRESS__"             23 20 30 0 \
     "use_spglib:"         24 1  "$__SCHA_USESPGLIB__"               24 20 30 0 \
    2>&1 1>&3)

# close fd
exec 3>&-
}

# ------------------------------------------------------------------------------
# relax parameters
__RELAX_TYPE__="relax"
__RELAX_NCONFIGS__=1000
__RELAX_MAX_POP__=1
__RELAX_START_POP__=1
__RELAX_SAVE_ENSEMBLE__="../ensemble_data_test"
__RELAX_GENERATE_FIRST_ENSEMBLE__=".false."
__RELAX_TARGET_PRESSURE__=0
__RELAX_FIXVOLUME__=".false."
__RELAX_BULK_MODULUS__=15

display_Relax_input() {
  # open fd
  exec 3>&1

  # Store data to $VALUES variable
  RELAX_VALUES=$(dialog --ok-label "Submit" \
	  --backtitle "The Stochastic Self-Consistent Harmonic Approximation (SSCHA)" \
	  --title "Relax parameters" \
	  --form "Type and n_confings are mandatory." \
    15 50 0 \
	   "type:"               1 1	"$__RELAX_TYPE__"      	             1 20 30 0 \
	   "n_configs:"          2 1	"$__RELAX_NCONFIGS__"  	             2 20 30 0 \
	   "max_pop_id:"         3 1	"$__RELAX_MAX_POP__"	               3 20 30 0 \
	   "start_pop:"          4 1	"$__RELAX_START_POP__"	             4 20 30 0 \
     "ensemble_datadir:"   5 1  "$__RELAX_SAVE_ENSEMBLE__"           5 20 30 0 \
     "generate_ensemble:"  6 1  "$__RELAX_GENERATE_FIRST_ENSEMBLE__" 6 20 30 0 \
     "target_pressure:"    7 1  "$__RELAX_TARGET_PRESSURE__"         7 20 30 0 \
     "fix_volume:"         8 1  "$__RELAX_FIXVOLUME__"               8 20 30 0 \
     "bulk_modulus:"       9 1  "$__RELAX_BULK_MODULUS__"            9 20 30 0 \
    2>&1 1>&3)

# close fd
exec 3>&-
}

# ------------------------------------------------------------------------------
# Calculator parameters

__KPTS_HEAD__=" 2 2 1 "
__KOFF_HEAD__=" 0 0 1 "
__DISABLE_CHECK__=.true.
__CALCULATOR_TYPE__="quantum-espresso"
__BINARY__=.false.
__QE_ALLOWED_KEYS__=


display_Calculator_input() {
  # open fd
  exec 3>&1

  # Store data to $VALUES variable
  CALCULATOR_VALUES=$(dialog --ok-label "Submit" \
	  --backtitle "The Stochastic Self-Consistent Harmonic Approximation (SSCHA)" \
	  --title "Calculator parameters" \
	  --form "program and k_points are mandatory." \
    15 50 0 \
	   "program:"            1 1	"$__CALCULATOR_TYPE__"               1 20 30 0 \
	   "k points:"           2 1	"$__KPTS_HEAD__"  	                 2 20 30 0 \
	   "k offset:"           3 1	"$__KOFF_HEAD__"	                   3 20 30 0 \
	   "binary:"             4 1	"$__BINARY__"	                       4 20 30 0 \
     "disable check:"      5 1  "$__DISABLE_CHECK__"                 5 20 30 0 \
     "QE extra keys:"      6 1  "$__QE_ALLOWED_KEYS__"               6 20 30 0 \
    2>&1 1>&3)

# close fd
exec 3>&-
}
# ------------------------------------------------------------------------------
# Utilities parameters

__UTILS_SAVEFREQ_FILENAME__="frequencies.dat"
__UTILS_SAVERHO_FILENAME__=
__UTILS_LOCKMODE_START__=
__UTILS_LOCKMODE_END__=
__UTILS_FREEMODE_START__=30
__UTILS_FREEMODE_END__=36
__UTILS_PROJECT_DYN__=.true.
__UTILS_PROJECT_STRUCTURE__=.false.



display_Utilities_input() {
  # open fd
  exec 3>&1

  # Store data to $VALUES variable
  UTILITIES_VALUES=$(dialog --ok-label "Submit" \
	  --backtitle "The Stochastic Self-Consistent Harmonic Approximation (SSCHA)" \
	  --title "Utilities parameters" \
	  --form "Optional" \
    15 50 0 \
	   "save freq filename:" 1 1	"$__UTILS_SAVEFREQ_FILENAME__"       1 20 30 0 \
	   "save_rho_filename:"  2 1	"$__UTILS_SAVERHO_FILENAME__"  	     2 20 30 0 \
	   "mu_lock_start:"      3 1	"$__UTILS_LOCKMODE_START__"	         3 20 30 0 \
	   "mu_lock_end:"        4 1	"$__UTILS_LOCKMODE_END__"	           4 20 30 0 \
     "mu_free_start:"      5 1  "$__UTILS_FREEMODE_START__"          5 20 30 0 \
     "mu_free_end:"        6 1  "$__UTILS_FREEMODE_END__"            6 20 30 0 \
     "project_dyn:"        7 1  "$__UTILS_PROJECT_DYN__"             7 20 30 0 \
     "project_structure:"  8 1  "$__UTILS_PROJECT_STRUCTURE__"       8 20 30 0 \
    2>&1 1>&3)

# close fd
exec 3>&-
}
# ------------------------------------------------------------------------------

display_result() {
  dialog --title "$1" \
    --no-collapse \
    --msgbox "$result" 0 0
}

display_help() {
  dialog --backtitle "The Stochastic Self-Consistent Harmonic Approximation (SSCHA)" \
  --title "SSCHA Help" --no-collapse --textbox sscha_help.txt 0 0
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
      display_Scha_input
      ;;
    2 )
      display_Calculator_input
      ;;
    3 )
        display_Relax_input
      ;;
    4)
        display_Utilities_input
      ;;
    5)
        display_help
      ;;
    6)
        result=$(echo "Scha:\n $SCHA_VALUES";uptime)
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
