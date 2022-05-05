#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH -N 1
#SBATCH --time=00:01:00
#SBATCH --job-name="template"

node=`hostname`
echo "******************"
echo "job run at node " $node
echo "******************"
echo ""
#copies the directory from where you submited the job to lscratch
cp -r $SLURM_SUBMIT_DIR  /lscratch/$USER/$SLURM_JOB_ID 
cd /lscratch/$USER/$SLURM_JOB_ID
rm slurm*.out

##########################################################################

echo "put your jobs here"

#########################################################################


cd $SLURM_SUBMIT_DIR
echo "Making backup..."
mkdir BACKUP

for file in *; do
	if [ $file != slurm*out ]; then
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
	if [ $file != BACKUP ] && [ $file != slurm*out ]; then
		rm -r BACKUP/$file
	fi
done
rmdir BACKUP



echo "DONE"
