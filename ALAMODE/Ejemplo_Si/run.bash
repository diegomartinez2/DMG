#!/bin/bash

cp harm1.lammps tmp.lammps
lmp_serial < in.sw > log.lammps
cp XFSET XFSET.harm1

for ((i=1;i<=81;i++))
do
    num=`echo $i | awk '{printf("%02d",$1)}'`
    cp cubic${num}.lammps tmp.lammps
    ./lmp_mpi_chimes < in.sw > log.lammps${num}
    cp XFSET XFSET.cubic${num}
done
