#!/bin/bash

np=40           #number of cpus
POPULATION=1    #population index
IONS=54         #number of atoms in the supercells
NCONFSSHA=300   #number of configurations in the sscha ensemble
#mpirun vasp_std > stdout
for i in `seq 1 $NCONFSSHA`; do
    echo `date` >> timing
    cp POSCAR_$i POSCAR
#    mpirun -np $np vasp_std > stdout
    ~/VASP/vasp.6.3.0/bin/vasp_std > stdout
    grep "energy  without entropy" OUTCAR  >> energies
    grep "forces" -A $IONS vasprun.xml > forces/forces_population$POPULATION_$i.dat
    rm POSCAR
    mv OUTCAR OUTCAR_$i
    echo `date` >> timing
done
#########################################################################
