#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  minimize.py
#
#  Copyright 2023 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
import sys, os, argparse #getopt

import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib import cm
import numpy as np
# ----------
# Funciones
# ----------
def sscha_run(POPULATION=1, N_RANDOM=100, SUPERCELL= (2,2,2), T=50, NQIRR=10):
    # Define input variables

    # NQIRR = 10                    # The number of irreducible q points in the q-point grid
    # T = 50                        # The temperature at which we want to perform the calculation in Kelvin
    # SUPERCELL = (2,2,2)           # The size of the supercell (or the q point grid)
    # N_RANDOM = 100                # The number of configurations that will be created
    # POPULATION = 1                # The population to generate

    # Load the dynamical matrices that generated the ensemble

    namefile='pop'+str(POPULATION)+'/dyn/dynq'
    if (POPULATION != 1):
        namefile='dyn_end_population'+str(POPULATION-1)+'_'
    dyn = CC.Phonons.Phonons(namefile, NQIRR)
    dyn.ForcePositiveDefinite()
    dyn.Symmetrize()
    # We make a copy of the starting dynamica matrices

    #dyn_0 = dyn

    # Prepare the stochastic weights for the SSCHA minimization

    folder_with_ensemble='pop'+str(POPULATION)+'/data'
    ensemble = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL)
    ensemble.load(folder_with_ensemble, population = POPULATION, N = N_RANDOM)
    #ensemble.update_weights(dyn, T)
    # Define the minimization

    minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

    # Define the steps for the centroids and the force constants

    #minimizer.min_step_dyn = 0.5 #0.005         # The minimization step on the dynamical matrix
    #minimizer.min_step_struc = 0.5 #0.05        # The minimization step on the structure
    #minimizer.kong_liu_ratio = 0.5         # The parameter that estimates whether the ensemble is still good
    #minimizer.meaningful_factor = 0.000001 # How much small the gradient should be before I stop?
    minimizer.set_minimization_step(0.01)
    minimizer.enforce_sum_rule=True
    #ensemble.update_weights

    # Let's start the minimization

    minimizer.init()

    # ioinfo=sscha.Utilities.IOInfo()
    # ioinfo.SetupSaving("minim_info")

    minimizer.run()

    # Let's make some plot on the evolution of the minimization

    name_for_data = 'population'+str(POPULATION)+'_minimization.dat'
    name_for_plot = 'population'+str(POPULATION)+'_minimization.pdf'

    minimizer.plot_results(save_filename = name_for_data, plot = False)

    fig, axs = plt.subplots(2,2)

    data = np.loadtxt(name_for_data)

    axs[0,0].plot(data[:,0],data[:,1])
    axs[0,1].plot(data[:,0],data[:,3], label='gradient')
    axs[0,1].plot(data[:,0],data[:,4], label='error')
    axs[1,0].plot(data[:,0],data[:,5], label='gradient')
    axs[1,0].plot(data[:,0],data[:,6], label='error')
    axs[1,1].plot(data[:,0],data[:,7])

    axs[0,0].set_xlabel('Step')
    axs[1,0].set_xlabel('Step')
    axs[0,1].set_xlabel('Step')
    axs[1,1].set_xlabel('Step')

    axs[0,0].set_ylabel(r'$F$')
    axs[0,1].set_ylabel(r'$\partial F/\partial \Phi$')
    axs[1,0].set_ylabel(r'$\partial F/\partial \mathcal{R}$')
    axs[1,1].set_ylabel(r'$N_{eff}$')

    axs[0,1].legend(frameon=False)
    axs[1,0].legend(frameon=False)

    plt.savefig(name_for_plot)

    # Let's print the final quantities
    # dyn_sc = minim.dyn.GenerateSupercellDyn( minim.ensemble.supercell )
    # # Diagonalize
    # w, pols = dyn_sc.DyagDinQ(0)
    # self.total_freqs.append(w)
    # np.savetxt(freq_name, self.total_freqs, header = "Time vs Frequencies")

    minimizer.finalize()

    # Let's print the final dynamical matrices

    namefile='dyn_end_population'+str(POPULATION)+'_'
    minimizer.dyn.save_qe(namefile)

    path = './pop'+str(POPULATION+1)+'/dyn'
    # Check whether the specified path exists or not
    isExist = os.path.exists(path)
    if not isExist:
       # Create a new directory because it does not exist
       os.makedirs(path)
    #os.makedirs('./pop'+str(POPULATION+1)+'/dyn')
    namefile='./pop'+str(POPULATION+1)+'/dyn/dynq'
    minimizer.dyn.save_qe(namefile)
    #save frequencies to file
    f = open("Frequencies"+str(POPULATION+1)+"_cm_1.txt", "a")
    w_s, pols = minimizer.dyn.DiagonalizeSupercell()
    f.write(" ".join(["{:.4f}".format(w * CC.Units.RY_TO_CM) for w in  w_s]))
    f.write("\n")
    f.close()
    #print convergence of calculation
    print ("SSCHA Converged? = {}".format(minimizer.is_converged()))
    return 0
    # harm_dyn = CC.Phonons.Phonons("harmonic_dyn", nqirr =3) # We load harmonic_dynX with X = 1,2,3
    # # Now we can print all the phonon frequencies
    # w_s, pols = harm_dyn.DiagonalizeSupercell()
    # print (" ".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in  w_s]))

def main(argv):
    # try:
    #     opts, args = getopt.getopt(argv,"hi:o:",["ipop=","iconf=","icell=","itemp=","inqirr="])
    # except getopt.GetoptError:
    #     print ('minimize.py -pop <population> -nconf <number of configurations> -cell <cell size> -temp <Temperature> -nqirr <number of irreducible q points>')
    # for opt, arg in opts:
    #     if opt == '-h':
    #         print ('minimnize.py -pop <population> -nconf <number of configurations> -cell <cell size)> -temp <Temperature> -nqirr <number of irreducible q points>')
    #         sys.exit()
    #     elif opt in ("-pop", "--ipop"):
    #         POPULATION = arg
    #     elif opt in ("-nconf", "--iconf"):
    #         N_RANDOM = arg
    #     elif opt in ("-cell", "--icell"):
    #         SUPERCELL_SIZE = arg
    #     elif opt in ("-temp", "--itemp"):
    #         T = arg
    #     elif opt in ("-nqirr", "--inqirr"):
    #         NQIRR = arg
    # SUPERCELL = (SUPERCELL_SIZE SUPERCELL_SIZE,SUPERCELL_SIZE,SUPERCELL_SIZE)
    argParser = argparse.ArgumentParser()
    argParser.add_argument("-pop", "--POPULATION", type=int, help="Population number")
    argParser.add_argument("-nconf", "--N_RANDOM", type=int, help="Number of configurations")
    argParser.add_argument("-cell", "--SUPERCELL", type=int, help="Cell size")
    argParser.add_argument("-temp", "--T", type=int, help="Temperature")
    argParser.add_argument("-nqirr", "--NQIRR", type=int, help="Number of irreducible q points")
    args = argParser.parse_args()
    SUPERCELL = (args.SUPERCELL,args.SUPERCELL,args.SUPERCELL)
    sscha_run(args.POPULATION, args.N_RANDOM, SUPERCELL, args.T, args.NQIRR)

if __name__ == "__main__":
   main(sys.argv[1:])
