#!/bin/bash
#SBATCH --job-name=VASP_300            # Job name
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=diego.martinez@ehu.eus  # Where to send mail
#!#SBATCH -p test                             # queue
#SBATCH --nodes=1                           # Run all processes on a single node
#SBATCH --ntasks=40                         # Number of processes
#SBATCH --time=44:30:00                     # Time limit hrs:min:sec
#SBATCH --output=job.log                    # Standard output and error log

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
############################################################################
module load VASP
#np=4           #number of cpus
POPULATION=$1    #population index
IONS=320         #number of atoms in the supercells
NCONFSSHA=512   #number of configurations in the sscha ensemble
#mpirun vasp_std > stdout
mkdir forces
for i in `seq 1 $NCONFSSHA`; do
    echo "--run--" $i
    echo `date` >> timing
    cp POSCAR_$i POSCAR
#    mpirun -np $np ~/VASP/vasp.6.3.0/bin/vasp_std > stdout
    mpirun -np 40 vasp_std > stdout
    grep "energy  without entropy" OUTCAR | awk '{print $NF}' >> energies
    #grep "energy  without entropy" OUTCAR  >> energies
    grep "forces" -A $IONS vasprun.xml > forces/forces_population$POPULATION'_'$i.dat
    rm POSCAR
    mv OUTCAR OUTCAR_$i
    echo `date` >> timing
done
#########################################################################
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
        if [ $file != BACKUP ] && [ $file != slurm*out ] && [ -e BACKUP/$file ]; then
                rm -r BACKUP/$file
        fi
done
rmdir BACKUP

echo "DONE"
