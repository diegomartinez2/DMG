Calculations of second-order phase transitions with the SSCHA
=============================================================

In this chapter we provide ready to use examples to calculate second-order phase transitions with SSCHA calculations.

Structural instability: calculation of the Hessian
--------------------------------------------------

This tutorial explains how to search for structural instabilities with a SSCHA calculation.

The SSCHA method provides a complete theoretical framework to study second-order phase transitions for structural instabilities.

According to Landau’s theory of second-order phase transitions, a phase transition occurs when the free energy curvature around the high-symmetry structure on the direction of the order parameter becomes negative:
[image]

For structural phase transitions, the order parameter is associated to phonon atomic displacements $$\frac{\partial^2 F}{\partial R_a \partial R_b}$$. So we just need to calculate the Free energy Hessian. the SSCHA provides an analytical equation for the free energy Hessian, derived by Raffaello Bianco in the work Bianco et. al. Phys. Rev. B 96, 014111 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.014111>.
The free energy curvature can be written in matrix form as:

$$\frac{\partial^2 F}{\partial {R_a}\partial {R_b}} = \Phi_{ab} + \sum_{cdef} \stackrel{(3)}{\Phi}_{acd}[1 - \Lambda\stackrel{(4)}{\Phi}]^{-1}_{cdef} \stackrel{(3)}{\Phi}_{efb}$$

Fortunately, this complex equation can be evaluated from the ensemble with a simple function call:

.. code-block:: python
   ensemble.get_free_energy_hessian()

Lets see a practical example, first we calculate the SSCHA dynamical matrix for the PbTe:

Lets do something different an write it as an object (Object Oriented Program).

.. code-block:: python
  class PbTe_initial(object):
    def __init__(self,fichero_ForceFields,fichero_dyn,nqirr,configuraciones,sobol,sobol_scatter):
        # Load the dynamical matrix for the force field
        self.ff_dyn = CC.Phonons.Phonons(fichero_ForceFields, 3)

        # Setup the forcefield with the correct parameters
        self.ff_calculator = ff.Calculator.ToyModelCalculator(self.ff_dyn)
        self.ff_calculator.type_cal = "pbtex"
        self.ff_calculator.p3 = 0.036475
        self.ff_calculator.p4 = -0.022
        self.ff_calculator.p4x = -0.014

        # Initialization of the SSCHA matrix
        self.dyn_sscha = CC.Phonons.Phonons(fichero_dyn, nqirr)
        # Flip the imaginary frequencies into real ones
        self.dyn_sscha.ForcePositiveDefinite()
        # Apply the ASR and the symmetry group
        self.dyn_sscha.Symmetrize()
        self.configuraciones=configuraciones
        self.sobol = sobol
        self.sobol_scatter = sobol_scatter

    def ensambla(self,T):
        self.ensemble = sscha.Ensemble.Ensemble(self.dyn_sscha, T0 = T, supercell = self.dyn_sscha.GetSupercell())
        # Detect space group
        symm=spglib.get_spacegroup(self.dyn_sscha.structure.get_ase_atoms(), 0.005)
        print('Initial SG = ', symm)


    def minimiza(self,fichero_frecuencias,fichero_matriz):
        self.minim = sscha.SchaMinimizer.SSCHA_Minimizer(self.ensemble)

        # Lets setup the minimization on the fourth root
#        self.minim.root_representation = "root4" # Other possibilities are 'normal' and 'sqrt'

        # To work correctly with the root4, we must deactivate the preconditioning on the dynamical matrix
#        self.minim.precond_dyn = False

        # Probemos con self.neglect_symmetries= true
#        self.minim.neglect_symmetries = True    # *test*

        # Now we setup the minimization parameters
        # Since we are quite far from the correct solution, we will use a small optimization step
        self.minim.min_step_dyn = 0.5 # If the minimization ends with few steps (less than 10), decrease it, if it takes too much, increase it

        # We decrease the Kong-Liu effective sample size below which the population is stopped
        self.minim.kong_liu_ratio = 0.5 # Default 0.5
        # Pedimos que minimize la estructura (eso espero) *test*
#        self.minim.minim_struct = True

#        self.relax = sscha.Relax.SSCHA(self.minim,
#                          ase_calculator = self.ff_calculator,
#                          N_configs = 100,
#                          max_pop = 200)
        self.relax = sscha.Relax.SSCHA(self.minim,
                          ase_calculator = self.ff_calculator,
                          N_configs = self.configuraciones,
                          max_pop = 50)

        # Setup the custom function to print the frequencies at each step of the minimization
        self.io_func = sscha.Utilities.IOInfo()
        self.io_func.SetupSaving(fichero_frecuencias) # The file that will contain the frequencies is frequencies.dat

        # Now tell relax to call the function to save the frequencies after each iteration
        # CFP stands for Custom Function Post (Post = after the minimization step)
        self.relax.setup_custom_functions(custom_function_post = self.io_func.CFP_SaveFrequencies)
        # Finalmente hacemos todos los calculos de busqueda de la energia libre.
        self.relax.relax(sobol = self.sobol, sobol_scramble = self.sobol_scatter)
#        self.relax.relax(sobol = False)
#        self.relax.vc_relax(static_bulk_modulus="recalc",restart_from_ens = True, fix_volume = True, stress_numerical = True)
        #self.relax.vc_relax(static_bulk_modulus=40, fix_volume = False)

        # Save the final dynamical matrix
        self.relax.minim.dyn.save_qe(fichero_matriz)
        # Detect space group
        symm=spglib.get_spacegroup(self.relax.minim.dyn.structure.get_ase_atoms(), 0.005)
        print('New SG = ', symm)
        view(self.relax.minim.dyn.structure.get_ase_atoms())

    def dibuja(self,fichero):
        # Setup the interactive plotting mode
        #plt.ion()

        # Lets plot the Free energy, gradient and the Kong-Liu effective sample size
        self.relax.minim.plot_results()

#        frequencies = np.loadtxt(fichero)
        frequencies = np.loadtxt("{}.freqs".format(fichero)) # se ha cambiado el formato?????
        N_steps, N_modes = frequencies.shape

        #For each frequency, we plot it [we convert from Ry to cm-1]
        plt.figure(dpi = 120)
        for i_mode in range(N_modes):
            plt.plot(frequencies[:, i_mode] * CC.Units.RY_TO_CM)
        plt.xlabel("Steps")
        plt.ylabel("Frequencies [cm-1]")
        plt.title("Evolution of the frequencies")
        plt.tight_layout()
        #plt.show()
        plt.savefig('Step_Freq.png')[...]

.. code-block:: python
  # Lets import all the sscha modules
  import cellconstructor as CC
  import cellconstructor.Phonons
  import sscha, sscha.Ensemble

  # We load the SSCHA dynamical matrix for the PbTe (the one after convergence)
  dyn_sscha = CC.Phonons.Phonons("dyn_sscha", nqirr = 3)

  # Now we load the ensemble
  ensemble = sscha.Ensemble.Ensemble(dyn_sscha, T0 = 1000, supercell=dyn_sscha.GetSupercell())
  ensemble.load("data_ensemble_final", N = 100, population = 5)

  # If the sscha matrix was not the one used to compute the ensemble
  # We must update the ensemble weights
  # We can also use this function to simulate a different temperature.
  ensemble.update_weights(dyn_sscha, T = 1000)

  # ----------- COMPUTE THE FREE ENERGY HESSIAN -----------
  dyn_hessian = ensemble.get_free_energy_hessian()
  # -------------------------------------------------------

  # We can save the free energy hessian as a dynamical matrix in quantum espresso format
  dyn_hessian.save_qe("free_energy_hessian")

We can then print the frequencies of the hessian. If an imaginary frequency is present, then the system wants to spontaneosly break the high symmetry phase.

The frequencies in the free energy hessian are temperature dependent.

*****

Plot the Hessian phonon dispersion
----------------------------------
