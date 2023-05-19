#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SSCHA_exercise.py
#
# Import the cellconstructor stuff
import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

# Import the modules of the force field
import fforces as ff
import fforces.Calculator

# Import the modules to run the sscha
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.Utilities

import spglib
from ase.visualize import view

# Import Matplotlib to plot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import timeit

class SnTe_initial(object):
  def __init__(self,file_ForceFields,file_dyn,nqirr,configurations,sobol,sobol_scatter):
      # Load the dynamical matrix for the force field
      self.ff_dyn = CC.Phonons.Phonons(file_ForceFields, 3)

      # Setup the forcefield with the correct parameters
      self.ff_calculator = ff.Calculator.ToyModelCalculator(self.ff_dyn)
      self.ff_calculator.type_cal = "pbtex"
      self.ff_calculator.p3 = 0.036475
      self.ff_calculator.p4 = -0.022
      self.ff_calculator.p4x = -0.014

      # Initialization of the SSCHA matrix
      self.dyn_sscha = CC.Phonons.Phonons(file_dyn, nqirr)
      # Flip the imaginary frequencies into real ones
      self.dyn_sscha.ForcePositiveDefinite()
      # Apply the ASR and the symmetry group
      self.dyn_sscha.Symmetrize()
      # Set some parameters and flags
      self.configurations=configurations
      self.sobol = sobol
      self.sobol_scatter = sobol_scatter

  def ensemble_sscha(self,T):
      self.ensemble = sscha.Ensemble.Ensemble(self.dyn_sscha, T0 = T, supercell = self.dyn_sscha.GetSupercell())
      # Detect space group
      symm=spglib.get_spacegroup(self.dyn_sscha.structure.get_ase_atoms(), 0.005)
      print('Initial SG = ', symm)


  def minimizing(self,file_frequencies,file_matrix):
      self.minim = sscha.SchaMinimizer.SSCHA_Minimizer(self.ensemble)

      # Lets setup the minimization on the fourth root
      #self.minim.root_representation = "root4" # Other possibilities are 'normal' and 'sqrt'

      # To work correctly with the root4, we must deactivate the preconditioning on the dynamical matrix
      #self.minim.precond_dyn = False

      # Now we setup the minimization parameters
      # Since we are quite far from the correct solution, we will use a small optimization step
      self.minim.min_step_dyn = 0.5 # If the minimization ends with few steps (less than 10), decrease it, if it takes too much, increase it

      # We decrease the Kong-Liu effective sample size below which the population is stopped
      self.minim.kong_liu_ratio = 0.5 # Default 0.5
      # We relax the structure
      self.relax = sscha.Relax.SSCHA(self.minim,
                        ase_calculator = self.ff_calculator,
                        N_configs = self.configurations,
                        max_pop = 50)

      # Setup the custom function to print the frequencies at each step of the minimization
      self.io_func = sscha.Utilities.IOInfo()
      self.io_func.SetupSaving(file_frequencies) # The file that will contain the frequencies is frequencies.dat

      # Now tell relax to call the function to save the frequencies after each iteration
      # CFP stands for Custom Function Post (Post = after the minimization step)
      #self.relax.setup_custom_functions(custom_function_post = self.io_func.CFP_SaveFrequencies)
      self.relax.setup_custom_functions(custom_function_post = self.io_func.CFP_SaveAll)
      # Finalmente hacemos todos los calculos de busqueda de la energia libre.
      self.relax.relax(sobol = self.sobol, sobol_scramble = self.sobol_scatter)
      #self.relax.vc_relax(static_bulk_modulus=40, fix_volume = False)

      # Save the final dynamical matrix
      self.relax.minim.dyn.save_qe(file_matrix)
      # Detect space group
      symm=spglib.get_spacegroup(self.relax.minim.dyn.structure.get_ase_atoms(), 0.005)
      print('New SG = ', symm)
      view(self.relax.minim.dyn.structure.get_ase_atoms())

class Search_instabilities(object):
    def __init__(self,files_ForceFields,files_dyn,nqirr):
        # Load the dynamical matrix for the force field
        self.ff_dyn = CC.Phonons.Phonons(files_ForceFields, 3)

        # Setup the forcefield with the correct parameters
        self.ff_calculator = ff.Calculator.ToyModelCalculator(self.ff_dyn)
        self.ff_calculator.type_cal = "pbtex"
        self.ff_calculator.p3 = 0.036475
        self.ff_calculator.p4 = -0.022
        self.ff_calculator.p4x = -0.014

        # Initialization of the SSCHA matrix
        self.dyn_sscha = CC.Phonons.Phonons(files_dyn, nqirr)
        self.dyn_sscha.ForcePositiveDefinite()

        # Apply also the ASR and the symmetry group
        self.dyn_sscha.Symmetrize()

    def load_dyn(self,File_final_dyn,nqirr):
        # The SSCHA dynamical matrix is needed (the one after convergence)
        # We reload the final result (no need to rerun the sscha minimization)
        self.dyn_sscha_final = CC.Phonons.Phonons(File_final_dyn, nqirr)

    def ensemble_sscha(self,T):
        # We reset the ensemble
        self.ensemble = sscha.Ensemble.Ensemble(self.dyn_sscha_final, T0 = T, supercell = self.dyn_sscha_final.GetSupercell())

        # We need a bigger ensemble to properly compute the hessian
        # Here we will use 10000 configurations
        self.ensemble.generate(5000, sobol = True, sobol_scramble = False)
        #self.ensemble.generate(10000, sobol = False)
        #We could also load the ensemble with ensemble.load("data_ensemble_final", N = 100, population = 5)

    def calculate(self):
        # We now compute forces and energies using the force field calculator
        self.ensemble.get_energy_forces(self.ff_calculator, compute_stress = False)

    def hessian(self,T):
        print("Updating the importance sampling...")
        # If the sscha matrix was not the one used to compute the ensemble
        # We must update the ensemble weights
        # We can also use this function to simulate a different temperature.
        self.ensemble.update_weights(self.dyn_sscha_final, T)
        # ----------- COMPUTE THE FREE ENERGY HESSIAN -----------
        print("Computing the free energy hessian...")
        self.dyn_hessian = self.ensemble.get_free_energy_hessian(include_v4 = False) # We neglect high-order four phonon scattering
        #self.dyn_hessian = self.ensemble.get_free_energy_hessian(include_v4 = True,
        #                                          get_full_hessian = True,verbose = True) # Full calculus
        # We can save the free energy hessian as a dynamical matrix in quantum espresso format
        self.dyn_hessian.save_qe("hessian")
        # -------------------------------------------------------
        # We calculate the frequencies of the hessian:
        w_hessian, pols_hessian = self.dyn_hessian.DiagonalizeSupercell()

        # Print all the frequency converting them into cm-1 (They are in Ry)
        print("\n".join(["{:16.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in w_hessian]))

class Hessian_Vs_Temperature(object):
    def __init__(self,T0,temperatures_i,files_ForceFields,configurations,sobol,sobol_scatter):
        # Load the dynamical matrix for the force field
        self.ff_dyn = CC.Phonons.Phonons(files_ForceFields, 3)

        # Setup the forcefield with the correct parameters
        self.ff_calculator = ff.Calculator.ToyModelCalculator(self.ff_dyn)
        self.ff_calculator.type_cal = "pbtex"
        self.ff_calculator.p3 = 0.036475
        self.ff_calculator.p4 = -0.022
        self.ff_calculator.p4x = -0.014
        # Define the temperatures, from 50 to 300 K, 6 temperatures
        #self.temperatures = np.linspace(50, 300, 6)
        self.temperatures = temperatures_i

        self.lowest_hessian_mode = []
        self.lowest_sscha_mode = []

        # Perform a simulation at each temperature
        self.t_old = T0
        self.configurations = configurations
        self.sobol = sobol
        self.sobol_scatter = sobol_scatter

    def cycle_T(self,Files_final_dyn,nqirr):
        for Temperature in self.temperatures:
            # Load the starting dynamical matrix
            self.dyn = CC.Phonons.Phonons(Files_final_dyn.format(int(self.t_old)), nqirr)

            # Prepare the ensemble
            self.ensemble = sscha.Ensemble.Ensemble(self.dyn, T0 = Temperature, supercell = self.dyn.GetSupercell())

            # Prepare the minimizer
            self.minim = sscha.SchaMinimizer.SSCHA_Minimizer(self.ensemble)
            self.minim.min_step_struc = 0.05
            self.minim.min_step_dyn = 0.002
            self.minim.kong_liu_ratio = 0.5
            self.minim.meaningful_factor = 0.000001
            #minim.root_representation = "root4"
            #minim.precond_dyn = False
            #self.minim.minim_struct = True
            #self.minim.neglect_symmetries = True
            self.minim.enforce_sum_rule = True  # Lorenzo's solution to the error

            # Prepare the relaxer (through many population)
            self.relax = sscha.Relax.SSCHA(self.minim, ase_calculator = self.ff_calculator, N_configs=self.configurations, max_pop=20)

            # Relax
            self.relax.relax(sobol = self.sobol, sobol_scramble = self.sobol_scatter)
            #self.relax.relax()

            # Save the dynamical matrix
            self.relax.minim.dyn.save_qe(Files_final_dyn.format(int(Temperature)))

            # Detect space group
            symm=spglib.get_spacegroup(self.relax.minim.dyn.structure.get_ase_atoms(), 0.005)
            print('Current SG = ', symm,' at T=',int(Temperature))

            # Recompute the ensemble for the hessian calculation
            self.ensemble = sscha.Ensemble.Ensemble(self.relax.minim.dyn, T0 = Temperature, supercell = self.dyn.GetSupercell())
            self.ensemble.generate(self.configurations, sobol = self.sobol, sobol_scramble = self.sobol_scatter)
            self.ensemble.get_energy_forces(self.ff_calculator, compute_stress = False) #gets the energies and forces from ff_calculator

            #update weights!!!
            self.ensemble.update_weights(self.relax.minim.dyn, Temperature)
            # Get the free energy hessian
            dyn_hessian = self.ensemble.get_free_energy_hessian(include_v4 = False) #free energy hessian as in Bianco paper 2017
            dyn_hessian.save_qe("hessian_T{}_".format(int(Temperature)))

            # Get the lowest frequencies for the sscha and the free energy hessian
            w_sscha, pols_sscha = self.relax.minim.dyn.DiagonalizeSupercell() #dynamical matrix
            # Get the structure in the supercell
            superstructure = self.relax.minim.dyn.structure.generate_supercell(self.relax.minim.dyn.GetSupercell()) #

            # Discard the acoustic modes
            acoustic_modes = CC.Methods.get_translations(pols_sscha, superstructure.get_masses_array())
            w_sscha = w_sscha[~acoustic_modes]

            self.lowest_sscha_mode.append(np.min(w_sscha) * CC.Units.RY_TO_CM) # Convert from Ry to cm-1

            w_hessian, pols_hessian = dyn_hessian.DiagonalizeSupercell() #recomputed dyn for hessian
            # Discard the acoustic modes
            acoustic_modes = CC.Methods.get_translations(pols_hessian, superstructure.get_masses_array())
            w_hessian = w_hessian[~acoustic_modes]
            self.lowest_hessian_mode.append(np.min(w_hessian) * CC.Units.RY_TO_CM) # Convert from Ry to cm-1
            #print ("\n".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in pols_hessian]))
            #exit()

            self.t_old = Temperature
        # We prepare now the file to save the results
        freq_data = np.zeros( (len(self.temperatures), 3))
        freq_data[:, 0] = self.temperatures
        freq_data[:, 1] = self.lowest_sscha_mode
        freq_data[:, 2] = self.lowest_hessian_mode

        # Save results on file
        np.savetxt("{}_hessian_vs_temperature.dat".format(self.configurations), freq_data, header = "T [K]; SSCHA mode [cm-1]; Free energy hessian [cm-1]")

    def draw_figure(self):
        hessian_data = np.loadtxt("{}_hessian_vs_temperature.dat".format(self.configurations))

        plt.figure(dpi = 120)
        plt.plot(hessian_data[:,0], hessian_data[:,1], label = "Min SCHA freq", marker = ">")
        plt.plot(hessian_data[:,0], hessian_data[:,2], label = "Free energy curvature", marker = "o")
        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("Temperature [K]")
        plt.ylabel("Frequency [cm-1]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}_Temp_Freq.png'.format(self.configurations))
        #plt.show()

        plt.figure(dpi = 120)
        plt.plot(hessian_data[:,0], np.sign(hessian_data[:,2]) * hessian_data[:,2]**2, label = "Free energy curvature", marker = "o")
        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("Temperature [K]")
        plt.ylabel("$\omega^2$ [cm-2]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}_Temp_Omeg.png'.format(self.configurations))
        #plt.show()


def main(args):
  print ("Running the SSCHA HANDS-ON-SESSION-3")
  #Setting the variables:
  #Setting the temperature in Kelvin:
  Temperature = 0
  #Setting the temparature range in Kelvin (6 steps from 50K to 300K):
  Temperature_i = np.linspace(50, 300, 6)
  #Setting the number of configurations:
  configuration_number = 50
  #Setting the names and location of the files:
  Files_ForceFields = "ffield_dynq"
  Files_dyn_SnTe = "ffield_dynq"
  #Set the number of irreducible q (reated to the supercell size):
  nqirr = 3
  #Setting the frequencies output file:
  File_frequencies = "frequencies.dat"
  #Setting the dynamical matrix output filename:
  File_final_dyn = "final_sscha_T{}_".format(int(Temperature))
  sobol = False
  sobol_scatter = False

  #We can comment this if we already made the sscha minimization
  #Calculus = SnTe_initial(Files_ForceFields,Files_dyn_SnTe,nqirr,configuration_number,sobol,sobol_scatter)
  #Calculus.ensemble_sscha(Temperature)
  #Calculus.minimizing(File_frequencies,File_final_dyn.format(int(Temperature)))

  #Now we can search for structural instabilities:
  #Unstable = Search_instabilities(Files_ForceFields,Files_dyn_SnTe,nqirr)
  #Unstable.load_dyn(File_final_dyn.format(Temperature),nqirr)
  #Unstable.ensemble_sscha(Temperature)
  #Unstable.calculate()
  #Unstable.hessian(Temperature)

  #
  HessianVsTemperature = Hessian_Vs_Temperature(Temperature,Temperature_i,Files_ForceFields,configuration_number,sobol,sobol_scatter)
  HessianVsTemperature.cycle_T(File_final_dyn,nqirr)
  HessianVsTemperature.draw_figure()

  return 0


if __name__ == '__main__':
  import sys
  sys.exit(main(sys.argv))
