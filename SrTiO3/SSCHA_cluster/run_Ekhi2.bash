#!/bin/bash
#SBATCH --job-name=tetragonal_2            # Job name
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=diego.martinez@ehu.eus  # Where to send mail
#SBATCH -p long                             # queue
#SBATCH --nodes=1                           # Run all processes on a single node
#SBATCH --ntasks=40                         # Number of processes
#SBATCH --time=16-00:00:00                     # Time limit hrs:min:sec
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
module load ASE/3.22.1-foss-2022a
source ~/SSCHA_1.4.1/sscha-1.4.1/bin/activate
export PYTHONPATH=~/.local/software/sscha_1.4.1/lib/python3.10/site-packages:$PYTHONPATH
export PATH=~/.local/software/sscha_1.4.1/bin:$PATH
echo "============================"
echo "Run SSCHA bash"
#np=4           #number of cpus
POPULATION=$1    #population index
IONS=320         #number of atoms in the supercells
NCONFSSHA=16384   #number of configurations in the sscha ensemble
#mpirun vasp_std > stdout
bash run_local.sh $POPULATION
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
