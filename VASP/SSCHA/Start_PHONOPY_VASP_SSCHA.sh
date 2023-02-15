#!/bin/bash
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
0.5 0.5 0.5" >> POSCAR
phonopy -d --dim='4 4 4'
cp POSCAR POSCAR_UNITCELL
cp SPOSCAR POSCAR
cp INCAR.harmonic INCAR
echo "VASP calculation"
phonopy --fc vasprun.xml
cp POSCAR_UNITCELL POSCAR
echo "interface"
cp INCAR.sc ./pop1/vasp/INCAR
cp POTCAR ./pop1/vasp/POTCAR
cp ML_FF ./pop1/vasp/ML_FF
nano run.bash
echo "Running VASP"
