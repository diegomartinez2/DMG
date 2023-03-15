#!/bin/bash

np=4           #number of cpus
POPULATION=1    #population index
IONS=40         #number of atoms in the supercells
NCONFSSHA=128   #number of configurations in the sscha ensemble
#mpirun vasp_std > stdout

for i in `seq 1 $NCONFSSHA`; do
    echo "--run--" $i
    echo `date` >> timing
    cp POSCAR_$i POSCAR
#    mpirun -np $np ~/VASP/vasp.6.3.0/bin/vasp_std > stdout
    ~/VASP/vasp.6.3.0/bin/vasp_std > stdout
    grep "energy  without entropy" OUTCAR | awk '{print $NF}' >> energies
    #grep "energy  without entropy" OUTCAR  >> energies
    grep "forces" -A $IONS vasprun.xml > forces/forces_population$POPULATION'_'$i.dat
    rm POSCAR
    mv OUTCAR OUTCAR_$i
    echo `date` >> timing
done
#########################################################################
