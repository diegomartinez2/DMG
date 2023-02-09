#!/bin/bash

#SBATCH --job-name=SrTiO3_VASP_4x           # Job name
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=diego.martinez@ehu.eus  # Where to send mail
#########SBATCH -p test                             # queue
#SBATCH --nodes=1                           # Run all processes on a single node
##SBATCH --ntasks=40                         # Number of processes
#SBATCH --time=1-23:00:00                   # Time limit hrs:min:sec
#SBATCH --output=job.log                    # Standard output and error log
#SBATCH --cpus-per-task=40

module load intel/2021a

node=`hostname`
echo "******************"
echo "job run at node " $node
echo "******************"
echo ""
#copies the directory from where you submited the job to lscratch

if [[ ! -e /lscratch/$USER ]]; then
    mkdir /lscratch/$USER
fi

cp -r $SLURM_SUBMIT_DIR  /lscratch/$USER/$SLURM_JOB_ID
cd /lscratch/$USER/$SLURM_JOB_ID
export NPROCS=$SLURM_NTASKS
rm slurm*.out

##########################################################################
module load VASP
np=40           #number of cpus
IONS=54         #number of atoms in the supercells
NCONFSSHA=300   #number of configurations in the sscha ensemble
#mpirun vasp_std > stdout
for i in {1..$NCONFSSHA}; do
    echo 'date' >> timing
    cp POSCAR_$i POSCAR
    mpirun -np $np vasp_std > stdout
    grep "energy  without entropy" OUTCAR  >> energies
    grep "forces" -A $IONS vasprun.xml > forces/forces_population1_$i.dat
    rm POSCAR
    mv OUTCAR OUTCAR_$i
    echo 'date' >> timing
done
#########################################################################


cd $SLURM_SUBMIT_DIR
echo "Making backup..."
mkdir BACKUP

for file in *; do
        if [ $file != BACKUP ] && [ $file != slurm*out ]; then
                mv $file BACKUP/$file
        fi
done

echo "Copying files from /lscratch..."
cp -r /lscratch/$USER/$SLURM_JOB_ID/* $SLURM_SUBMIT_DIR/

echo "Deleting files from /lscratch..."
rm -r /lscratch/$USER/$SLURM_JOB_ID

echo "Deleting the BACKUP..."
cd $SLURM_SUBMIT_DIR
for file in *; do
        if [$file != BACKUP] && [$file != slurm*out] && [-e BACKUP/$file]; then
                rm -r BACKUP/$file
        fi
done
rmdir BACKUP

echo "DONE"
