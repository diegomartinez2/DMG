#!/bin/bash
if [[-z $1]]; then
  echo "Use Start_PHONOPY_VASP_SSCHA 0 'N N N' for the harmonic calculation of NxNxN supercell (write it like '2 2 2' with the tildes).
  Use Start_PHONOPY_VASP_SSCHA 1 for the first VASP force/energy calculation."
elif  [[$1 -eq 0]]; then
  echo "SrTiO3
 3.9048
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
 Sr O Ti
 1 3 1
direct
0.0 0.0 0.0
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0
0.5 0.5 0.5" > POSCAR
  echo " IBRION = 6
 ML_FF_LMLFF = .TRUE.
 ML_FF_ISTART = 2" > INCAR.harmonic
 echo " ML_FF_LMLFF = .TRUE.
 ML_FF_ISTART = 2" > INCAR.sc
#phonopy -d --dim='4 4 4'
  phonopy -d --dim=$2
  cp POSCAR POSCAR_UNITCELL
  cp SPOSCAR POSCAR
  cp INCAR.harmonic INCAR
  echo "VASP harmonic calculation"; sleep 300
elif [[$1 -eq 1]]; then
  phonopy --fc vasprun.xml
  cp POSCAR_UNITCELL POSCAR
  echo "interface"
  cp INCAR.sc ./pop1/vasp/INCAR
  cp POTCAR ./pop1/vasp/POTCAR
  cp ML_FF ./pop1/vasp/ML_FF
  nano run.bash
  echo "Running VASP"
else
  echo "I do not understand you."
fi
