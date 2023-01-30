#!/bin/bash
#
#   USAGE: script to extract energy from each folder.

#
s=0
echo -e "************************************************************************* \e[92m"
echo "*************************************************************************"
echo -e "*************************************************************************\e[0m"
echo -e "\e[0m" | awk '{printf "%15s %15s %15s\n", "Directory", "ENERGY (eV)", "VOLUME (A**3)"}'
for d in */;
do
  E = 'grep "free  energy ML TOTEN  =" OUTCAR | tail -1 | awk '{printf "%f", $6 }''
#  echo "$E" | awk '{printf "%6.6f \n ", $1* 0.0367493 } ' >> .tmp1  #Bohr units
  echo "$E" | awk '{printf "%6.6f \n ", $1* 0.0734986176 } ' >> .tmp1
  paste .tmp1 > energy

  echo -e " \e[0m"
  unset s
