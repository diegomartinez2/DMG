Apendix: the Ekhi cluster
=========================

For this SSCHA School we will provide access to the local CFM cluster named Ekhi.

ekhi.cfm.ehu.es
---------------

Ekhi cluster, designed specifically for novel Quantum ESPRESSO calculations, is composed of 28 computing nodes with two Xeon Cascade Lake-SP 6230 processors (40 computing cores) and 96 GB of memory in each node, with an Infiniband FDR interconnection network, giving a total of 1120 cores and 2.7 TB of memory.

You can view the cluster information with *sinfo*:

.. code-block:: bash

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

.. code-block:: bash

  #!/bin/bash
  #SBATCH --job-name=Test_ESPRESSO            # Job name
  #SBATCH --mail-type=NONE                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
  #SBATCH --mail-user=                        # Where to send mail
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

.. code-block:: bash

  sbatch run.sh


A typical usage of this cluster from a SSCHA script code includes:

.. code-block:: python

  #-----------------------------------------------------------------------
  username = user_name   # Put here your login name for the cluster.
  pseudo = {"Sr": "Sr.pbesol-spn-kjpaw_psl.1.0.0.UPF",
            "Ti": "Ti.pbesol-spn-kjpaw_psl.1.0.0.UPF",
            "O" : "O.pbesol-n-kjpaw_psl.1.0.0.UPF"}
  input_params = {"tstress" : True, # Print the stress in the output
          "tprnfor" : True, # Print the forces in the output
          "tstress" : True, #output stresses
          "ecutwfc" : 70,  #The wavefunction energy cutoff for plane-waves (Ry)
          "ecutrho" : 700, # The density energy cutoff (Ry)
          "mixing_beta" : 0.4,  # The mixing parameter in the self-consistent calculation
          "conv_thr" : 1e-9,    # The energy convergence threshold (Ry)
          "degauss" : 0.03,  # Smearing temperature (Ry)
  #                "smearing" : "mp",
          "pseudo_dir" : "./pseudo/",
          "occupations" : "fixed", #smearing or fixed (fixed for insulators with a gap; gaussian smearing for metals; )
          "disk_io" : "none"}

  k_points = (8,8,8) # The k points grid (you can alternatively specify a kspacing)
  k_offset = (1,1,1) # The offset of the grid (can increase convergence)

  self.espresso_calc = Espresso(pseudopotentials = pseudo, input_data = input_params,
                  kpts = k_points, koffset = k_offset)
  my_hpc = sscha.Cluster.Cluster(pwd = None)
  # We setup the connection info
  my_hpc.hostname = "{}@ekhi.cfm.ehu.es".format(username) # The command to connect via ssh to the cluster (pippo@login.cineca.marconi.it)
  my_hpc.workdir = "/scratch/{}/my_calculation".format(username) # the directory in which the calculations are performed

  # Now we need to setup the espresso
  # First we must tell the cluster where to find him:
  my_hpc.binary = "pw.x -npool NPOOL -i  PREFIX.pwi > PREFIX.pwo"
  # Then we need to specify if some modules must be loaded in the submission script
  my_hpc.load_modules = """
  # Here this is a bash script at the beginning of the submission
  # We can load modules

  module load QuantumESPRESSO
  export OMP_NUM_THREADS=1
  """

  # All these information are independent from the calculation
  # Now we need some more specific info, like the number of processors, pools and other stuff
  my_hpc.n_cpu = 40 # We will use 32 processors
  my_hpc.n_nodes = 1 #In 1 node
  my_hpc.n_pool = 10 # This is an espresso specific tool, the parallel CPU are divided in 4 pools

  # We can also choose in how many batch of jobs we want to submit simultaneously, and how many configurations for each job
  my_hpc.batch_size = 10
  my_hpc.job_number = 10
  # In this way we submit 10 jobs, each one with 10 configurations (overall 100 configuration at time)

  # We give 25 seconds of timeout
  my_hpc.set_timeout(25)

  # We can specify the time limit for each job,
  my_hpc.time = "03:00:00" # 5 minutes

  # Create the working directory if none on the cluster
  # And check the connection
  my_hpc.setup_workdir()
  #-----------------------------------------------------------------------

Then we can use in relax with:

.. code-block:: python

  relax = sscha.Relax.SSCHA(minim, ase_calculator = espresso_calc, N_configs=configurations, max_pop=20, cluster = my_hpc)
