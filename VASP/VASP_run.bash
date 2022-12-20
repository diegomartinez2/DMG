#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=500:00:00
#PBS -N Si
#PBS -q batch
#PBS -j oe
cd $PBS_O_WORKDIR
#Linux PART
ulimit -s unlimited
ulimit -n unlimited
ulimit -l unlimited
EXEC=/home/shuang/softwares/VASP/bin/vasp_std
NCORE=`cat $PBS_NODEFILE | wc -l`
mpirun -np $NCORE $EXEC >log 2>error
