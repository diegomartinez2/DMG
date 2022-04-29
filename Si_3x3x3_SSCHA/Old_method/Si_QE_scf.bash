#!/bin/bash
#
#SBATCH --job-name=Si_QE_scf_mpi
#SBATCH --output=res_mpi.txt
#
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=diego.martinez@ehu.eus # Where to send mail	
#
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=4		    # 6 cpus (best ratios 4,5,8,10)
#SBATCH --time=1-23:59:00		    # exageramos el tiempo (ya lo ajustamos luego)
#SBATCH --mem-per-cpu=1gb	    # ya veremos cuanto se necesita
#SBATCH --output=array_%A-%a.out    # Standard output and error log
#SBATCH --array=1-50                # Array range ???

module load QuantumESPRESSO

FILES=(./scf*.in)

srun mpirun -np $SLURM_NTASKS pw.x < ${FILES[$SLURM_ARRAY_TASK_ID]} > ${FILES[$SLURM_ARRAY_TASK_ID]}.out
