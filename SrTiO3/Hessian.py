import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble

# We load the SSCHA dynamical matrix for the PbTe (the one after convergence)
#dyn_sscha = CC.Phonons.Phonons("dyn_sscha", nqirr = 4)
dyn_sscha = CC.Phonons.Phonons("sscha_T300_dyn", nqirr = 4)
# Now we load the ensemble
ensemble = sscha.Ensemble.Ensemble(dyn_sscha, T0 = 300, supercell=dyn_sscha.GetSupercell())
#ensemble.load("data_ensemble_final", N = 1024, population = 5)
ensemble.generate(1000,sobol = True)
# If the SSCHA matrix was not the one used to compute the ensemble
# We must update the ensemble weights
# We can also use this function to simulate a different temperature.
ensemble.update_weights(dyn_sscha, 1000)

# ----------- COMPUTE THE FREE ENERGY HESSIAN -----------
dyn_hessian = ensemble.get_free_energy_hessian()
# -------------------------------------------------------

# We can save the free energy hessian as a dynamical matrix in quantum espresso format
dyn_hessian.save_qe("free_energy_hessian")
