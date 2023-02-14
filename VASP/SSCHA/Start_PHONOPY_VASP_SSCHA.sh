#!/bin/bash
echo "SrTiO3 Unit cell
3.9048 # cell parameter in angstrom
Sr 0 0 0
O 0.0 0.5 0.5
O 0.5 0.0 0.5
O 0.5 0.5 0.0
Ti 0.5 0.5 0.5"
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
