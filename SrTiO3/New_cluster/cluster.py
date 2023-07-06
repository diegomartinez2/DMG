import cellconstructor as CC, cellconstructor.Phonons
import sscha
import sscha.Cluster
import sys, os
def configure_cluster(cluster_workdir = "SrTiO3_workdir"):
    cluster = sscha.Cluster.Cluster(hostname = "user11@ekhi.cfm.ehu.es")
