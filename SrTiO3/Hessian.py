import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble

# We load the SSCHA dynamical matrix for the PbTe (the one after convergence)
#dyn_sscha = CC.Phonons.Phonons("dyn_sscha", nqirr = 4)
dyn_sscha = CC.Phonons.Phonons("sscha_T300_dyn", nqirr = 4)
dyn_sscha2 = CC.Phonons.Phonons("dyn_pop4_", nqirr = 4) # <- la de la ultima poblacion
# Now we load the ensemble
ensemble = sscha.Ensemble.Ensemble(dyn_sscha2, T0 = 300, supercell=dyn_sscha.GetSupercell()) #<- la de la ultima poblacion??? 
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

w_hessian, pols_hessian = dyn_hessian.DiagonalizeSupercell()
    # Discard the acoustic modes
ss = dyn_hessian.structure.generate_supercell(dyn_hessian.GetSupercell())

acoustic_modes = CC.Methods.get_translations(pols_hessian, ss.get_masses_array())

w_hessian = w_hessian[~acoustic_modes]
##lowest_hessian_mode.append(np.min(w_hessian) * CC.Units.RY_TO_CM) # Convert from Ry to cm-1
# Print all the frequency converting them into cm-1 (They are in Ry)
print("\n".join(["{:16.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in w_hessian]))
