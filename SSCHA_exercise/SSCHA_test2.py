#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SSCHA_test2.py
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
# ---------------------------
# Importación de los módulos
# ---------------------------
from __future__ import print_function
from __future__ import division

import numpy
import matplotlib
from matplotlib import pylab, mlab, pyplot
np = numpy
plt = pyplot

from IPython.core.pylabtools import figsize, getfigs

from pylab import *
from numpy import *


import sys,os
#importamos el metodo de calculo de las fuerzas Born-Oppenhaimer (Quantum-espresso)
import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view
#importamos los modulos base del SSCHA (cellconstructor)
import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons

#importamos el motor SSCHA
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax

import spglib


# -----------
# Constantes
# -----------

# ------------------------------
# Clases y Funciones utilizadas
# ------------------------------
class Mi_sscha(object):
    def __init__(self):
        self.pseudo = {"H": "H.pbe-rrkjus_psl.1.0.0.UPF",
                    "La" : "La.pbe-spfn-rrkjus_psl.1.0.0.UPF"}
        self.input_params = {"tstress" : True, # Print the stress in the output
                    "tprnfor" : True, # Print the forces in the output
                    "ecutwfc" : 35,  #The wavefunction energy cutoff for plane-waves (Ry)
                    "ecutrho" : 350, # The density energy cutoff (Ry)
                    "mixing_beta" : 0.2,  # The mixing parameter in the self-consistent calculation
                    "conv_thr" : 1e-9,    # The energy convergence threshold (Ry)
                    "degauss" : 0.02,  # Smearing temperature (Ry)
                    "smearing" : "mp",
                    "pseudo_dir" : ".",
                    "occupations" : "smearing",
                    "disk_io" : "none"}

        self.k_points = (8,8,8) # The k points grid (you can alternatively specify a kspacing)
        self.k_offset = (1,1,1) # The offset of the grid (can increase convergence)

        self.espresso_calc = Espresso(pseudopotentials = pseudo, input_data = input_params, kpts = k_points, koffset = k_offset)

    def uno(self):
        # We now load the dynamical matrix
        self.dyn = CC.Phonons.Phonons("dyn")
        self.dyn.Symmetrize() #Enforce the sum rule
        # We prepare the ensemble
        self.ensemble = sscha.Ensemble.Ensemble(dyn, T0 = 0, supercell = dyn.GetSupercell())
        # We prepare the sscha minimizer
        self.minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
        # We set up the minimization parameters
        self.minim.min_step_dyn = 0.05     # The minimization step on the dynamical matrix
        self.minim.min_step_struc = 0.05   # The minimization step on the structure
        self.minim.kong_liu_ratio = 0.5     # The parameter that estimates whether the ensemble is still good
        self.minim.gradi_op = "all" # Check the stopping condition on both gradients
        self.minim.meaningful_factor = 0.2 # How much small the gradient should be before I stop?

    def cluster(self):
        # Here we prepare a cluster
        # Here we configure the cluster object MARCONI
        import sscha.Cluster
        self.my_hpc = sscha.Cluster.Cluster(pwd = None)
        # We setup the connection info
        self.my_hpc.hostname = "ekhi" # The command to connect via ssh to the cluster
        #my_hpc.account_name = "IscrB_COMRED" # The name of the project for the computation
        self.my_hpc.workdir = "/scratch/lorenzo/my_calculation" # the directory in which the calculations are performed
        # Now we need to setup the espresso
        # First we must tell the cluster where to find him:
        self.my_hpc.binary = "pw.x -npool NPOOL -i  PREFIX.pwi > PREFIX.pwo"
        # Then we need to specify if some modules must be loaded in the submission script
        self.my_hpc.load_modules = """
# Here this is a bash script at the beginning of the submission
# We can load modules

module load QuantumESPRESSO
export OMP_NUM_THREADS=1
"""

        # All these information are independent from the calculation
        # Now we need some more specific info, like the number of processors, pools and other stuff
        self.my_hpc.n_cpu = 32 # We will use 32 processors
        self.my_hpc.n_nodes = 1 #In 1 node
        self.my_hpc.n_pool = 16 # This is an espresso specific tool, the parallel CPU are divided in 4 pools

        # We can also choose in how many batch of jobs we want to submit simultaneously, and how many configurations for each job
        self.my_hpc.batch_size = 20
        self.my_hpc.job_number = 20
        # In this way we submit 10 jobs, each one with 10 configurations (overall 100 configuration at time)

        # We give 25 seconds of timeout
        self.my_hpc.set_timeout(25)

        # We can specify the time limit for each job,
        self.my_hpc.time = "00:10:00" # 5 minutes

        # Create the working directory if none on the cluster
        # And check the connection
        self.my_hpc.setup_workdir()
    def relaja(self):
        # Decomment the following line if you did not set up the cluster
        #my_hpc = None

        self.relax = sscha.Relax.SSCHA(minim, ase_calculator = self.espresso_calc,
            N_configs = 400,
            max_pop = 20,
            save_ensemble = True,
            cluster = self.my_hpc)

        print ("The original spacegroup is:", spglib.get_spacegroup(self.dyn.structure.get_ase_atoms(), 0.05))

        self.space_groups = []
        self.relax.setup_custom_functions(custom_function_post = self.print_spacegroup)

        # Now we can run the calculation!!!
        # In this case we fix the volume (we optimize lattice parameters)
        # But you can also fixe the target pressure (as done in the commented line)
        if not os.path.exists("ensembles"):
            os.mkdir("ensembles")
        self.relax.vc_relax(fix_volume = True, static_bulk_modulus = 120, ensemble_loc = "ensembles")
        #self.relax.vc_relax(target_press = 120, static_bulk_modulus = 200, ensemble_loc = "ensembles")
        self.relax.minim.finalize()
        self.relax.minim.plot_results()

        spglib.get_spacegroup(self.relax.minim.dyn.structure.get_ase_atoms(), 0.05)
        spglib.get_spacegroup(self.dyn.structure.get_ase_atoms(), 0.1)
        view(self.relax.minim.dyn.structure.get_ase_atoms())
        self.relax.minim.dyn.structure.unit_cell


    def print_spacegroup(self):
        spgroup = spglib.get_spacegroup(self.minim.dyn.structure.get_ase_atoms(), 0.05)
        self.space_groups.append(spgroup)

        # We can save them in the output at each minimization step
        f = open("space_group.dat", "w")
        f.writelines(["{}) {}\n".format(i+1, x) for i,x in enumerate(self.space_groups)])
        f.close()

def main(args):
    calculo = Mi_sscha()
    calculo.uno()
    calculo.cluster()
    calculo.relaja()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
