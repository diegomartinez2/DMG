Ekhi
====

Ekhi cluster, designed specifically for novel Quantum ESPRESSO calculations, is composed of 28 computing nodes with two Xeon Cascade Lake-SP 6230 processors (40 computing cores) and 96 GB of memory in each node, with an Infiniband FDR interconnection network, giving a total of 1120 cores and 2.7 TB of memory.

You can view the cluster information with *sinfo*:

.. code:: bash

  PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
  all*         up 2-00:00:00      1    mix ekhi31
  all*         up 2-00:00:00     20  alloc ekhi[1-5,7-9,11-12,15-17,20-21,24-25,27-29]
  all*         up 2-00:00:00     10   idle ekhi[6,10,13-14,18-19,22-23,26,30]
  fat          up 14-00:00:0      1  alloc ekhi29
  fat          up 14-00:00:0      1   idle ekhi30
  long         up   infinite      4  alloc ekhi[1-4]
  test         up      30:00      1    mix ekhi31
  test         up      30:00     20  alloc ekhi[1-5,7-9,11-12,15-17,20-21,24-25,27-29]
  test         up      30:00     10   idle ekhi[6,10,13-14,18-19,22-23,26,30]


This is an example of a batch input for the cluster Ekhi:

.. code:: bash

  #!/bin/bash
  #SBATCH --job-name=Test_Espresso            # Job name
  #SBATCH --mail-type=ALL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
  #SBATCH --mail-user=diego.martinez@ehu.eus  # Where to send mail
  #SBATCH -p test                             # queue (in this example the queue for test)
  #SBATCH --nodes=1                           # Run all processes on a single node
  #SBATCH --ntasks=40                         # Number of processes
  #SBATCH --time=00:30:00                     # Time limit hrs:min:sec (for the test queue time max. = 30min.)
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


This batch is run with:

.. code:: bash

  sbatch run.sh
