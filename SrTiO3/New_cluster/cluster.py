import cellconstructor as CC, cellconstructor.Phonons
import sscha
import sscha.Cluster
import sys, os
def configure_cluster(cluster_workdir = "SrTiO3_workdir"):
    cluster = sscha.Cluster.Cluster(hostname = "user11@ekhi.cfm.ehu.es", pwd="RsAUUeVo")

    cluster.n_nodes = 1 #specifies the number of nodes
    cluster.custom_param["ntastks"] = 40

    cluster.time = "00:30:00" #total time
    cluster.n_cpu = 40 #specify how many processors to call quantum Espresso with
    cluster.n_pool = 4 #is the number of pools for the quantum espresso parallelization; it should be the greatest common divisor between the number of CPUs and K points.
    cluster.job_number = 10 #how many jobs will be submitted simultaneously (executed in parallel, but with queue time)
    cluster.batch_size = 1 #how many pw.x calculations to group in the same job (executed one after the other without queue time)

    scratch_workdir = os.path.join("/scratch/$USER/", cluster_workdir)
    cluster.workdir = scratch_workdir
    cluster.add_set_minus_x = True # Avoid the set -x
    cluster.load_modules = f"""

module load QuantumESPRESSO

export OMP_NUM_THREADS=1
"""

    cluster.setup_workdir()

    def cp_files(lbls):
        extrain = f"cd {scratch_workdir}\n"
        extraout = "sleep 1\n"
        for lbl in lbls:
            extrain += f"cp {home_workdir}/{lbl}.pwi {scratch_workdir}/\n"
            extraout += f"mv {scratch_workdir}/{lbl}.pwo {home_workdir}/\n"

        return extrain, extraout

    # Add the possibility to copy the input files
    #cluster.additional_script_parameters = cp_files
    #Check the communication
    if not cluster.CheckCommunication():
        raise ValueError("Impossible to connect to the cluster")

    return cluster
