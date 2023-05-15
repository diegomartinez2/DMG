==================
Hands-on-session 3
==================

Calculations of second-order phase transitions with the SSCHA
=============================================================

In this hands-on we provide ready to use examples to calculate second-order phase transitions with SSCHA calculations.

Structural instability: calculation of the Hessian
--------------------------------------------------

This tutorial explains how to search for structural instabilities with a SSCHA calculation.

The SSCHA method provides a complete theoretical framework to study second-order phase transitions for structural instabilities.

According to Landau’s theory of second-order phase transitions, a phase transition occurs when the free energy curvature around the high-symmetry structure on the direction of the order parameter becomes negative:

.. _fig-second_order:

.. figure:: figures/second_order.png
   :width: 400
   :alt: Second order.

   Landau’s theory of second-order phase transitions.

For structural phase transitions, the order parameter is associated to phonon atomic displacements:

.. math:: \frac{\partial^2 F}{\partial R_a \partial R_b}

So we just need to calculate the Free energy Hessian. the SSCHA provides an analytical equation for the free energy Hessian, derived by Raffaello Bianco in the work `Bianco et. al. Phys. Rev. B 96, 014111 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111>`_.
The free energy curvature can be written in matrix form as:

.. math:: \frac{\partial^2 F}{\partial {R_a}\partial {R_b}} = \Phi_{ab} + \sum_{cdef} \stackrel{(3)}{\Phi}_{acd}[1 - \Lambda\stackrel{(4)}{\Phi}]^{-1}_{cdef} \stackrel{(3)}{\Phi}_{efb}

Fortunately, this complex equation can be evaluated from the ensemble with a simple function call:

.. code-block:: python

   ensemble.get_free_energy_hessian()

Lets see a practical example, first we calculate the SSCHA dynamical matrix for the SnTe:

To speedup the calculations, we will use a force-field that can mimic the physics of ferroelectric transitions in FCC lattices.

Lets do something different an write it as an object (Object Oriented Program).

We begin from the bottom with the old python trick for code re-usability:

.. code-block:: python

  def main(args):
  #main code goes here
  return 0

  if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))

and now we put the top part of the code, we import some libraries:

.. code-block:: python

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

Now we need to calculate the SSCHA dynamical matrix. For that we can use this object:

.. code-block:: python

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
        # Finally we do all the free energy calculations.
        self.relax.relax(sobol = self.sobol, sobol_scramble = self.sobol_scatter)
        #self.relax.vc_relax(static_bulk_modulus=40, fix_volume = False)

        # Save the final dynamical matrix
        self.relax.minim.dyn.save_qe(file_matrix)
        # Detect space group
        symm=spglib.get_spacegroup(self.relax.minim.dyn.structure.get_ase_atoms(), 0.005)
        print('New SG = ', symm)
        view(self.relax.minim.dyn.structure.get_ase_atoms())


We've seen most of this before, so let's review what's there in detail:

1. First this object is initialized in "__init__" where the toy model potential is set for the next calculations.
   This force field needs the harmonic dynamical matrix to be initialized, and the higher order parameters.
   Finally, the dynamical matrix for the minimization is loaded and readied. Since we are studying a system that has a spontaneous symmetry breaking at low temperature, the harmonic dynamical matrices will have imaginary phonons. We must enforce phonons to be positive definite to start a SSCHA minimization.

2. The next function of this object "ensemble_sscha" just creates the ensembles for the specified temperature. As an extra, we also look for the space group of the structure.

3. Next comes the minimization function. In the "minimizing" function we can set the fourth root minimization, in which, instead of optimizing the auxiliary dynamical matrices themselves, we will optimize their fourth root.

   .. math:: \Phi = \left({\sqrt[4]{\Phi}}\right)^4

   This constrains the dynamical matrix to be positive definite during the minimization.
   Next the automatic relaxation is set with the option here to use the Sobol sequence for the ensemble generation.

   We also set a custom function to save the frequencies at each iteration, to see how they evolves. This is very useful to understand if the algorithm is converged or not.

   Then the dynamical matrix of the converged minimization is saved in a file, and finally we take a look at the space group and the structure.


Now we fill the main function with this new object:

.. code-block:: python

  def main(args):
    #Setting the variables:
    #Setting the temperature in Kelvin:
    Temperature = 0
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

    Calculus = SnTe_initial(Files_ForceFields,Files_dyn_SnTe,nqirr,configuration_number,sobol,sobol_scatter)
    Calculus.ensemble_sscha(Temperature)
    Calculus.minimizing(File_frequencies,File_final_dyn.format(int(Temperature)))

    return 0

This code will calculate the SSCHA minimization with the "ff_calculator". We cat use "sscha-plot-data.py" to take a look at the minimization.

.. code-block:: bash

    python sscha-plot-data.py frequencies.dat


Note: this force field model is not able to compute stress, as it is defined only at fixed volume, so we cannot use it for a variable cell relaxation.

**Now we can search for instabilities.**

If we have a very small mode in the SSCHA frequencies, it means that associated to that mode we have huge fluctuations. This can indicate an instability. However, to test this we need to compute the free energy curvature along this mode. This can be obtained in one shot thanks to the theory developed in `Bianco et. al. Phys. Rev. B 96, 014111. <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111>`_

For that we create another object, "Search_instabilities" to do the job.

[...]


.. code-block:: python

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

Lets see what this code do:

1. In the initialization function the "ff_calculator" toy potential is defined as we have seen in the previous object.

2. In "load_dyn" function, it will load the dynamical matrix calculated previously with the "ff_calculator" toy potential, so there is no need to calculate it again.

3. Then, as the Hessian calculation is more sensible, we generate a new ensemble with more configurations in "ensemble_sscha".
   To compute the hessian we will use an ensemble of 10000 configurations.
   Note here that we can use less if we use Sobol sequence or we can load a previously generated ensemble.

4. We now compute forces and energies using the force field calculator.

5. Finally the free energy hessian is calculated in the "hessian" function.
   We can choose if we neglect or not in the calculation the four phonon scattering process. Four phonon scattering processes require a huge memory allocation for big systems, that scales as (3⋅N)^4 with N the number of atoms in the supercell. Moreover, it may require also more configurations to converge.

   In almost all the systems we studied up to now, we found this four phonon scattering at high order to be negligible. We remark, that the SSCHA minimization already includes four phonon scattering at the lowest order perturbation theory, thus neglecting this term only affects combinations of one or more four phonon scattering with two three phonon scatterings (high order diagrams). For more details, see `Bianco et. al. Phys. Rev. B 96, 014111. <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111>`_

   We can then print the frequencies of the hessian. If an imaginary frequency is present, then the system wants to spontaneously break the high symmetry phase.

The frequencies in the free energy hessian are temperature dependent.

Lets put this object into the main function and calculate:

.. code-block:: python

  def main(args):
    #Setting the variables:
    #Setting the temperature in Kelvin:
    Temperature = 0
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
    Unstable = Search_instabilities(Files_ForceFields,Files_dyn_SnTe,nqirr)
    Unstable.load_dyn(File_final_dyn.format(Temperature),nqirr)
    Unstable.ensemble_sscha(Temperature)
    Unstable.calculate()
    Unstable.hessian(Temperature)

    return 0

We can look at the eigenmodes of the free energy hessian to check if we have imaginary phonons. If there are negative frequencies then we found an instability. You can check what happens if you include the fourth order.

Phase transition:
-----------------

Up to now we studied the system at T=0K and we found that there is an instability. However, we can repeat the minimization at many temperatures, and track the phonon frequency to see which is the temperature at which the system becomes stable.

.. code-block:: python

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

Lets see the code:

1. The initialization is similar to the one we did for the "Search_instabilities".

2. In "cycle_T" we condense in one function the calculation of the hessians in a loop for different temperatures. In the end, it searches for the lowest non acoustic frequency to save with the correspondent auxiliar sscha frequency.

3. Finally, in the last function of this object, we make a graphic output of the data.

Lets put this object into the main function and calculate:

.. code-block:: python

            def main(args):
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

We will simulate the temperatures up to room temperature (300 K) with steps of 50 K. Note, this will perform all the steps above 6 times, so it may take some minutes, depending on the PC (on a i3 from 2015, with one core, it took 2 hours).
If it takes too long you can reduce the number of steps by changing 'Temperature_i = np.linspace(50, 300, 6)'.

We are using only 50 configurations in the ensemble. Note that this makes a fast calculation but is a low number. How the calculation changes with the number of configurations?

As exercise, you can modify this "Hessian_Vs_Temperature" object by calling the "Search_instabilities" into the "cycle_T" function.



*****

Plot the Hessian phonon dispersion
----------------------------------
