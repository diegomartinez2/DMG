#!/bin/bash
#SBATCH --job-name=Spectral_function            # Job name
#SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=diego.martinez@ehu.eus  # Where to send mail
# #SBATCH --partition=fat                     # queue fat, long
#SBATCH --nodes=1                           # Run all processes on a single node
#SBATCH --ntasks=20                        # Number of processes
#SBATCH --time=2-00:00:00                   # Time limit hrs:min:sec
#SBATCH --output=job.log                    # Standard output and error log
# #SBATCH --men=272.1GB
# #SBATCH --men=755GB

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
source ~/SSCHA/sscha-foss-new-2023-11-29/bin/activate
module load ASE/3.22.1-foss-2022a
export PYTHONPATH=~/.local/software/sscha-2023-11-29/lib/python3.10/site-packages:$PYTHONPATH
export PATH=~/.local/software/sscha-2023-11-29/bin:$PATH
echo "============================"
echo "Run mpirun Hessian.py"
###mpirun python Hessian.py
#time python Hessian.py 2>&1 | tee Hessian_v3.txt
#time python Hessian_v4.py 2>&1 | tee Hessian_v4.txt
time python Spectral_function_calc.py 2>&1 | tee Spectral_function_calc.txt
echo "============================"
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
