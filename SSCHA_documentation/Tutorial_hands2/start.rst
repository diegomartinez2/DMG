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

.. code-blocl:: python
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
