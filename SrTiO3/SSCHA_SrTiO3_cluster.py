#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SSCHA_SrTiO3_cluster.py
#
#  Copyright 2022 Diego <diego@u038025>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# Import the sscha code
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax, sscha.Utilities

# Import the cellconstructor library to manage phonons
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.Structure, cellconstructor.calculators

# Import the force field of Gold
import ase, ase.calculators
from ase.calculators.emt import EMT

# Import numerical and general pourpouse libraries
import numpy as np, matplotlib.pyplot as plt
import sys, os


import cellconstructor.ForceTensor
import ase.dft.kpoints

import matplotlib.cm as cm
import matplotlib.colors as colors

import scipy, scipy.optimize

import sscha.Cluster

class Busca_inestabilidades(object):
    def __init__(self,fichero_ForceFields,fichero_dyn,nqirr):

        # Initialization of the SSCHA matrix
        self.dyn_sscha = CC.Phonons.Phonons(fichero_dyn, nqirr)
        self.dyn_sscha.ForcePositiveDefinite()

        # Apply also the ASR and the symmetry group
        self.dyn_sscha.Symmetrize()
    def load_dyn(self,Fichero_final_matriz_dinamica,nqirr):
        # We reload the final result (no need to rerun the sscha minimization)
        self.dyn_sscha_final = CC.Phonons.Phonons(Fichero_final_matriz_dinamica, nqirr)
    def ensambla(self,T):
        # We reset the ensemble
        self.ensemble = sscha.Ensemble.Ensemble(self.dyn_sscha_final, T0 = T, supercell = self.dyn_sscha_final.GetSupercell())

        # We need a bigger ensemble to properly compute the hessian
        # Here we will use 10000 configurations
        self.ensemble.generate(1024, sobol = True, sobol_scramble = False)
#        self.ensemble.generate(50, sobol = False)
#        self.ensemble.generate(1000,sobol = True)
    def calcula1(self):
        # We now compute forces and energies using the force field calculator
        self.ensemble.get_energy_forces(self.calculator, compute_stress = False) #test compute_stress = True no puede con este potencial...
    def hessiano(self,T):
        #self.dyn_hessian = self.ensemble.get_free_energy_hessian(include_v4 = False) # We neglect high-order four phonon scattering

        print("Updating the importance sampling...")
        self.ensemble.update_weights(self.dyn_sscha_final, T)

        print("Computing the free energy hessian...")
       #self.dyn_hessian = self.ensemble.get_free_energy_hessian(include_v4 = False) # We neglect high-order four phonon scattering
        self.dyn_hessian = self.ensemble.get_free_energy_hessian(include_v4 = True,
                                                  get_full_hessian = True,verbose = True) # Full calculus
        # We can save it
        self.dyn_hessian.save_qe("hessian")

        w_hessian, pols_hessian = self.dyn_hessian.DiagonalizeSupercell()

        # Print all the frequency converting them into cm-1 (They are in Ry)
        print("\n".join(["{:16.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in w_hessian]))

class Send_to_cluster(object):
    def __init__(self,hostname = 'diegom@ekhi.cfm.ehu.es', pwd = None,
           label = 'SrTiO3_', account_name = '', n_nodes = 1, n_cpu = 40,
           time = '00:20:00', n_pool = 20, workdir = '/scratch/diegom/SrTiO3',
           mpi_cmd=r"srun --mpi=pmi2 -n NPROC"):    #note the r"..." means that "..." is literate.
        self.cluster = sscha.Cluster.Cluster(hostname = hostname, pwd = pwd,  # Put the password in pwd if needed
            mpi_cmd = mpi_cmd)

        # Configure the submission strategy
        if (account_name != ''):
            self.cluster.account_name = account_name  # Name of the account on which to subtract nodes
        if (account_name == ''):
            self.cluster.use_account = False
        self.cluster.n_nodes = n_nodes            # Number of nodes requested for each job
        self.cluster.time = time                  # Total time requested for each job
        self.cluster.n_pool = n_pool              # Number of pools for the Quantum ESPRESSO calculation
        self.cluster.label = label                # label, change it if you are going to submit
                                                  # two different calculations in the same working directory
        # Here some custom parameters for the clusters
        # These are specific for daint, but you can easily figure out those for your machine
        #self.cluster.custom_params["--constraint"] = "gpu"      # Run on the GPU partition
        #self.cluster.custom_params["--ntasks-per-node"] = '2'
        #self.cluster.custom_params["--cpus-per-task"] = '6'

        # Since daint specify the partition with a custom option,
        # Lets remove the specific partition option of SLURM
        # Neither we want to specify the total number of cpus (automatically determined by the node)
        self.cluster.use_partition = False
        self.cluster.partition_name = ""
        self.cluster.use_cpu = True #False
        self.cluster.n_cpu = n_cpu


        # Now, we need to tell daint which modules to load to run quantum espresso
        # Also this is cluster specific, but very simple to figure it out for you
        self.cluster.load_modules = """
# Load the quantum espresso modules
## module load daint-gpu
module load QuantumESPRESSO

# Configure the environmental variables of the job
##export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NO_STOP_MESSAGE=1
## export CRAY_CUDA_MPS=1

## ulimit -s unlimited
        """

        # Now, what is the command to run quantum espresso on the cluster?
        self.cluster.binary = "pw.x -npool NPOOL -i PREFIX.pwi > PREFIX.pwo"  #<--- No need for mpirun it's set with mpi_cmd
        # NOTE that NPOOL will be replaced automatically with the cluster.n_pool variable


        # Let us setup the working directory (directory in which the jobs runs)
        #self.cluster.workdir =  "$SCRATCH/diegom/Gold_NVT_300k"       # <--- If this fails just use "/scratch/[...]" instead of $SCRATCH/[...]
        self.cluster.workdir = workdir
        self.cluster.setup_workdir()  # Login to the cluster and creates the working directory if it does not exist


        # Last but not least:
        #  How many jobs do you want to submit simultaneously?
        self.cluster.batch_size = 10

        #  Now many DFT calculation do you want to run inside each job?
        self.cluster.job_number = 5


class  SrTiO3_free_energy_ab_initio(object):
    def __init__(self):
        # Initialize the DFT (Quantum Espresso) calculator for SrTiO3
        # The input data is a dictionary that encodes the pw.x input file namelist
        input_data = {
            'control' : {
                # Avoid writing wavefunctions on the disk
                'disk_io' : 'None',
                # Where to find the pseudopotential
                #'pseudo_dir' : '.'
                'tstress' : True,
                'tprnfor' : True
            },
            'system' : {
                # Specify the basis set cutoffs
                'ecutwfc' : 60,   # Cutoff for wavefunction from 70
                'ecutrho' : 60*10, # Cutoff for the density ecutwfc*10
                # Information about smearing (it is a metal)
                'occupations' : 'fixed', # 'fixed' or 'smearing', smearing for conductors
                #'smearing' : 'mv',
                'degauss' : 0.03
            },
            'electrons' : {
                'conv_thr' : 1e-8,
                'conv_thr' : 1.e-09,
                'electron_maxstep' : 80,
                'mixing_beta' : 4.e-01
            }
        }

        # the pseudopotential for each chemical element
        # In this case O, Sr and Ti
        pseudopotentials = {'O' : 'O.pbesol-n-kjpaw_psl.1.0.0.UPF',
                            'Sr' : 'Sr.pbesol-spn-kjpaw_psl.1.0.0.UPF',
                            'Ti' : 'Ti.pbesol-spn-kjpaw_psl.1.0.0.UPF'}

        # the kpoints mesh and the offset
        kpts = (4,4,4)    #With the supercell the number of k points in effect are multiplied
        koffset = (0,0,0)

        # Specify the command to call quantum espresso
        command = 'mpirun -np 8 pw.x -i PREFIX.pwi > PREFIX.pwo'


        # Prepare the quantum espresso calculator
        self.calculator = CC.calculators.Espresso(input_data,
                                             pseudopotentials,
                                             command = command,
                                             kpts = kpts,
                                             koffset = koffset)

    def relax(self):

        TEMPERATURE = 300
        N_CONFIGS = 32
        MAX_ITERATIONS = 20
        START_DYN = 'harmonic_dyn'
        NQIRR = 10

        # Let us load the starting dynamical matrix
        SrTiO3_dyn = CC.Phonons.Phonons(START_DYN, NQIRR)
        SrTiO3_dyn.ForcePositiveDefinite()
        SrTiO3_dyn.Symmetrize()
        # Initialize the random ionic ensemble
        ensemble = sscha.Ensemble.Ensemble(SrTiO3_dyn, TEMPERATURE)

        # Initialize the free energy minimizer
        minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
        minim.set_minimization_step(0.01)

        # Initialize the NVT simulation
        mi_cluster = Send_to_cluster(hostname = 'diegom@ekhi.cfm.ehu.es', label = 'SrTiO3_', n_pool = 20, # n_pool must be a divisor of the k_points example 5x5x5=125 n_pool= 5 or 25
           n_cpu = 40, time = '00:20:00', mpi_cmd = 'mpirun -np NPROC' ) #test with 5 pools for QE; note reducing the n_pool reduces the memory usage in the k points calculation.
        relax = sscha.Relax.SSCHA(minim, self.calculator, N_configs = N_CONFIGS,
                        max_pop = MAX_ITERATIONS, cluster = mi_cluster.cluster,
                        save_ensemble = True)

        # Define the I/O operations
        # To save info about the free energy minimization after each step
        ioinfo = sscha.Utilities.IOInfo()
        ioinfo.SetupSaving("minim_info")
        relax.setup_custom_functions(custom_function_post = ioinfo.CFP_SaveAll)


        # Run the NVT simulation (save the stress to compute the pressure)
        relax.relax(get_stress = False, sobol = True, sobol_scatter = 0.0)

        # If instead you want to run a NPT simulation, use
        # The target pressure is given in GPa.
        #relax.vc_relax(target_press = 0)

        # You can also run a mixed simulation (NVT) but with variable lattice parameters
        #relax.vc_relax(fix_volume = True)

        # Now we can save the final dynamical matrix
        # And print in stdout the info about the minimization
        relax.minim.finalize()
        relax.minim.dyn.save_qe("sscha_T{}_dyn".format(TEMPERATURE))

def plot_dispersion_SrTiO3(PATH = "GX"):
    NQIRR = 10
    #CMAP = "Spectral_r"
    #PATH = "GX"
    #PATH = "GM"
    #PATH = "GR"

    N_POINTS = 1000

    SPECIAL_POINTS = {"G": [0,0,0],
                      "X": [0, 0, .5],
                      "M": [0, .5, .5],
                      "R": [.5, .5, .5]}

    # Load the harmonic and sscha phonons
    harmonic_dyn = CC.Phonons.Phonons('harmonic_dyn', NQIRR)
    sscha_dyn = CC.Phonons.Phonons('sscha_T300_dyn', NQIRR)

    # Get the band path
    qpath, data = CC.Methods.get_bandpath(harmonic_dyn.structure.unit_cell,
                                          PATH,
                                          SPECIAL_POINTS,
                                          N_POINTS)
    xaxis, xticks, xlabels = data # Info to plot correclty the x axis

    # Get the phonon dispersion along the path
    harmonic_dispersion = CC.ForceTensor.get_phonons_in_qpath(harmonic_dyn, qpath)
    sscha_dispersion = CC.ForceTensor.get_phonons_in_qpath(sscha_dyn, qpath)

    nmodes = harmonic_dyn.structure.N_atoms * 3

    # Plot the two dispersions
    plt.figure(dpi = 150)
    ax = plt.gca()

    for i in range(nmodes):
        lbl=None
        lblsscha = None
        if i == 0:
            lbl = 'Harmonic'
            lblsscha = 'SSCHA'

        ax.plot(xaxis, harmonic_dispersion[:,i], color = 'k', ls = 'dashed', label = lbl)
        ax.plot(xaxis, sscha_dispersion[:,i], color = 'r', label = lblsscha)

    # Plot vertical lines for each high symmetry points
    for x in xticks:
        ax.axvline(x, 0, 1, color = "k", lw = 0.4)
    ax.axhline(0, 0, 1, color = 'k', ls = ':', lw = 0.4)

    ax.legend()

    # Set the x labels to the high symmetry points
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

    ax.set_xlabel("Q path")
    ax.set_ylabel("Phonons [cm-1]")

    plt.tight_layout()
    plt.savefig("dispersion{}.png".format(PATH))
    plt.show()
    return 0

def main(args):
    SrTiO3_calculation = SrTiO3_free_energy_ab_initio()
    SrTiO3_calculation.relax()
    plot_dispersion_SrTiO3(PATH = "GX")
    plot_dispersion_SrTiO3(PATH = "GM")
    plot_dispersion_SrTiO3(PATH = "GR")
    #return 0
    #raise SystemExit
    #sys.exit()

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
