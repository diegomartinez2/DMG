#!/bin/bash
#SBATCH --job-name=Harmonic_VASP            # Job name
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=diego.martinez@ehu.eus  # Where to send mail
#SBATCH --partition=fat                             # queue
#SBATCH --nodes=1                           # Run all processes on a single node
##SBATCH --ntasks=40                         # Number of processes
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
mpirun vasp_std > stdout
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
