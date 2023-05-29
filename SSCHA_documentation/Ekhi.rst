Ekhi
====

Ekhi cluster, designed specifically for novel Quantum ESPRESSO calculations, is composed of 28 computing nodes with two Xeon Cascade Lake-SP 6230 processors (40 computing cores) and 96 GB of memory in each node, with an Infiniband FDR interconnection network, giving a total of 1120 cores and 2.7 TB of memory.

.. code:: bash

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
  module load QuantumESPRESSO
  ############################################################################
  echo "put your jobs here"
  mpirun -np 40 pw.x -npool 20 -i input_espresso.pwi > output_espresso.pwo
  ############################################################################


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
