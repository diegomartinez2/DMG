#!/bin/bash

#cp harm1.lammps tmp.lammps
#lmp_serial < in.sw > log.lammps
#cp XFSET XFSET.harm1

for ((i=1;i<=720;i++))
do
    num=`echo $i | awk '{printf("%03d",$1)}'`
    cp harm${num}.lammps tmp.lammps
    ./lmp_mpi_chimes < in.sw
    cp XFSET XFSET.harm${num}
done
for ((i=1;i<=16884;i++))
do
    num=`echo $i | awk '{printf("%05d",$1)}'`
    cp cubic${num}.lammps tmp.lammps
    ./lmp_mpi_chimes < in.sw
    cp XFSET XFSET.cubic${num}
done
