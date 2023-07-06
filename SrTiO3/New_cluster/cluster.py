import cellconstructor as CC, cellconstructor.Phonons
import sscha
import sscha.Cluster
import sys, os
def configure_cluster(cluster_workdir = "SrTiO3_workdir"):
    cluster = sscha.Cluster.Cluster(hostname = "user11@ekhi.cfm.ehu.es", pwd="RsAUUeVo")

    cluster.n_nodes = 1
    cluster.custom_param["ntastks"] = 40

    cluster.time = "00:30:00"
    cluster.n_cpu = 40
    cluster.n_pool = 4
    cluster.job_number = 10
    cluster.batch_size = 1

    scratch_workdir = os.path.join("/scratch/$USER/", cluster_workdir)
    cluster.workdir = scratch_workdir
    cluster.add_set_minus_x = True # Avoid the set -x
    cluster.load_modules = f"""

module load QuantumESPRESSO

export OMP_NUM_THREADS=1
"""

    cluster.setup_workdir()

    #Check the communication
    if not cluster.CheckCommunication():
        raise ValueError("Impossible to connect to the cluster")

    return cluster
