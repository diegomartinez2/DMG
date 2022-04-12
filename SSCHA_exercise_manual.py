#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
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
import sys,os

import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble
import sscha.SchaMinimizer, sscha.Relax

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib import cm
import numpy as np
# -------
# Clases
# -------
class Calculo(object):
    def __init__(self,Temperatura,NQirr,configuraciones):
        # Define input variables

        self.NQIRR = NQirr #3                     # The number of irreducible q points in the q-point grid
        self.T = Temperatura #200                       # The temperature at which we want to perform the calculation in Kelvin
        self.SUPERCELL = (2,2,2)           # The size of the supercell (or the q point grid)
        self.N_RANDOM = configuraciones #50                # The number of configurations that will be created
        self.POPULATION = 1                # The population to generate

        # define calculator
        self.ff_calculator = ff.Calculator.ToyModelCalculator(self.ff_dyn)
        self.ff_calculator.type_cal = "pbtex"
        self.ff_calculator.p3 = 0.036475
        self.ff_calculator.p4 = -0.022
        self.ff_calculator.p4x = -0.014

    def create_configurations(self):
        # Load the starting dynamical matrices in the 2x2x2 grid
        if self.POPULATION != 1:
            namefile='dyn_end_population'+str(self.POPULATION-1)+'_'  # We will start from the output dynamical matrices at the previous step
        else:
            namefile='harmonic_calculation/harmonic_dyn'  # We will start from the harmonic dynamical matrices at the first step


        self.dyn = CC.Phonons.Phonons(namefile, self.NQIRR)

        # Apply the sum rule and symmetries

        self.dyn.Symmetrize()

        # Flip the imaginary frequencies into real ones if there are imaginary phonon frequencies

        self.dyn.ForcePositiveDefinite()

        # We can print the frequencies to show the magic:
        #w_s, pols = dyn.DiagonalizeSupercell()
        #print ("\n".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in  w_s]))

        # Generate the ensemble

        self.ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL)
        self.ens.generate(N_RANDOM)

        # Save the ensemble

        namefile='population'+str(self.POPULATION)+'_ensemble'
        self.ens.save(namefile, self.POPULATION)

        # Prepare qe input files

        for i in range(self.N_RANDOM):
            bash_command = 'cat head_scf.in population'+str(self.POPULATION)+'_ensemble/scf_population'+str(self.POPULATION)+'_'+str(i+1)+'.dat > population'+str(self.POPULATION)+'_ensemble/scf_'+str(i+1)+'.in'
            os.system(bash_command)
    def calcula_qe(self):
        self.ens.compute_ensemble(ff_calculator)
        self.ens.save('data_enseble_ff',self.POPULATION)
    def extrae_energias(self):
        return 0
    def minimiza(self):
        # Load the dynamical matrices that generated the ensemble

        #namefile='dyn_start_population'+str(self.POPULATION)+'_'
        #self.dyn = CC.Phonons.Phonons(namefile, self.NQIRR)

        # We make a copy of the starting dynamica matrices

        dyn_0 = self.dyn

        # Prepare the stochastic weights for the SSCHA minimization

        folder_with_ensemble='population'+str(POPULATION)+'_ensemble'
        ensemble = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL)
        ensemble.load(folder_with_ensemble, population = POPULATION, N = N_RANDOM)

        # Define the minimization

        minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

        # Define the steps for the centroids and the force constants

        minimizer.min_step_dyn = 0.005         # The minimization step on the dynamical matrix
        minimizer.min_step_struc = 0.05        # The minimization step on the structure
        minimizer.kong_liu_ratio = 0.5         # The parameter that estimates whether the ensemble is still good
        minimizer.meaningful_factor = 0.000001 # How much small the gradient should be before I stop?

        # Let's start the minimization

        minimizer.init()
        minimizer.run()

        # Let's make some plot on the evolution of the minimization

        name_for_data = 'population'+str(POPULATION)+'_minimization.dat'
        #name_for_plot = 'population'+str(POPULATION)+'_minimization.pdf'

        minimizer.plot_results(save_filename = name_for_data, plot = False)

        #fig, axs = plt.subplots(2,2)

        #data = np.loadtxt(name_for_data)

        # axs[0,0].plot(data[:,0],data[:,1])
        # axs[0,1].plot(data[:,0],data[:,3], label='gradient')
        # axs[0,1].plot(data[:,0],data[:,4], label='error')
        # axs[1,0].plot(data[:,0],data[:,5], label='gradient')
        # axs[1,0].plot(data[:,0],data[:,6], label='error')
        # axs[1,1].plot(data[:,0],data[:,7])
        #
        # axs[0,0].set_xlabel('Step')
        # axs[1,0].set_xlabel('Step')
        # axs[0,1].set_xlabel('Step')
        # axs[1,1].set_xlabel('Step')
        #
        # axs[0,0].set_ylabel(r'$F$')
        # axs[0,1].set_ylabel(r'$\partial F/\partial \Phi$')
        # axs[1,0].set_ylabel(r'$\partial F/\partial \mathcal{R}$')
        # axs[1,1].set_ylabel(r'$N_{eff}$')
        #
        # axs[0,1].legend(frameon=False)
        # axs[1,0].legend(frameon=False)
        #
        # plt.savefig(name_for_plot)

        # Let's print the final quantities

        minimizer.finalize()

        # Let's print the final dynamical matrices

        namefile='dyn_end_population'+str(POPULATION)+'_'
        minimizer.dyn.save_qe(namefile)
    def Chequeo(self):
        self.running = not self.minim.is_converged() #para hacer "while running:" con paso a paso
        slef.POPULATION += 1
# ----------
# Funciones
# ----------

def main(args):
    temperatura = 0
    nquirr = 3
    numero_de_configuraciones = 10

    calcula = Calculo(temperatura,nquirr,numero_de_configuraciones)
    calcula.create_configurations()
    calcula.calcula_qe()
    calcula.minimiza()
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
