#Pre-processing crear la matriz dinámica
1. Crea el objeto con la estructura:

   structure = CC.Structure.Structure()

   1. Lee el archivo con la estructura en formato cif:

   structure.read_generic_file("Au.cif")

   2. Prepara el calculador:

   calculator = ase.calculators.emt()  #por ejemplo

   3. Relax structure with the calculator:

   relax = CC.calculators.Relax(structure, calculator)
   relaxed = relax.static_relax()

   4. Crea la matriz dinámica:

   dyn = CC.Phonons.compute_phonons_finite_displacements(structure, calculator, supercell = (4,4,4))

   5. Guarda la matriz dinámica:

   dyn.Symmetrize()
   dyn.save_qe(namefile)

# pasos minimizacion

1. Lee las matrices dinamicas:

    dyn = CC.Phonons.Phonons(namefile, NQIRR)

   1. Aplica la regla de la suma y simetriza:

    dyn.Symmetrize()

   2. Elimina las frecuencias imaginarias (forzando a que sean positivas definidas):

    dyn.ForcePositiveDefinite()

   * Como extra podemos ver las frecuencias_

   w_s, pols = dyn.DiagonalizeSupercell()

2. Genera el ensablado:

   1. ensemble = sscha.Ensemble.Ensemble(dyn, 1000, supercell = dyn.GetSupercell())

   2. ensemble.generate(N)

---
Ejecuta el calculo de las energias, fuerzas y el stress con un programa externo
---
Se repiten los pasos 1 y 2, pero cargando los resultados del calculo anterior

3. Define la minimizacion:

   minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

   1. Parametros de la minimizacion:

   minimizer.min_step_dyn = 0.005         # The minimization step on the dynamical matrix
   minimizer.min_step_struc = 0.05        # The minimization step on the structure
   minim.gradi_op = "all"                 # Check the stopping condition on both gradients
   minimizer.kong_liu_ratio = 0.5         # The parameter that estimates whether the ensemble is still good
   minimizer.meaningful_factor = 0.000001 # How much small the gradient should be before I stop?

Ahora hay dos opciones:
4. A: Paso a paso

   minimizer.init()
   minimizer.run()
   ...
   minimizer.finalize()

4. B: Relajación automatica:

   relax = sscha.Relax.SSCHA(minimizer,
                          ase_calculator = ff_calculator,
                          N_configs = 10000,
                          max_pop = 20)
   relax.relax()

## post processing: hesiano de Energia libre
%necesitamos actualizar los pesos!!!-> ens.update_weights(final_dyn, T)
* Matriz dinamica original:

   dyn = CC.Phonons.Phonons(DYN_PREFIX, NQIRR)

* Matriz dinamica actual:

   final_dyn = CC.Phonons.Phonons(FINAL_DYN, NQIRR)

* reensmblado:

   ens = sscha.Ensemble.Ensemble(dyn, Tg, dyn.GetSupercell())
   ens.load(DATA_DIR, POPULATION, N_RANDOM)

* Importante reajustar los pesos:

   ens.update_weights(final_dyn, T)

   dyn_hessian = ens.get_free_energy_hessian(include_v4 = INCLUDE_V4,
                                          get_full_hessian = True,
                                          verbose = True)
   dyn_hessian.save_qe(SAVE_PREFIX)

## post processing: spectral_function

#! Initialize the tensor3 object
#! We need 2nd FCs of the used grid to configure the supercell.
#! For example, we can use the sscha final auxiliary dynamical matrices
dyn = CC.Phonons.Phonons("dyn_end_population3_",3)
supercell = dyn.GetSupercell()
tensor3 = CC.ForceTensor.Tensor3(dyn.structure,
                                dyn.structure.generate_supercell(supercell),
                                supercell)

#! Assign the tensor3 values
d3 = np.load("d3_realspace_sym.npy")*2.0 # The 2 factor is because of units, needs to be passed to Ry
tensor3.SetupFromTensor(d3)

#! Center and apply ASR, which is needed to interpolate the third order force constant
tensor3.Center()
tensor3.Apply_ASR()

#! Print the tensor if you want, uncommenting the next line
#tensor3.WriteOnFile(fname="FC3",file_format='D3Q')

#! Calculate the spectral function at Gamma in the no-mode mixing approximation
#! keeping the energy dependence on the self-energy.
#
#! An interpolation grid needs to be used (and one needs to check convergence with
#! respect to it considering different values of the smearing)


#! interpolation grid
k_grid=[20,20,20]    

#
G=[0.0,0.0,0.0]

CC.Spectral.get_diag_dynamic_correction_along_path(dyn=dyn,
                                                   tensor3=tensor3,  
                                                   k_grid=k_grid,
                                                   q_path=G,
                                                   T = 200,                             # The temperature for the calculation  
                                                   e1=145, de=0.1, e0=0,                # The energy grid in cm-1
                                                   sm1=1.0, nsm=1, sm0=1.0,             # The smearing \eta for the analytic continuation
                                                   filename_sp = 'nomm_spectral_func')  # Output file name

#! Now perform the calculation of the spectral function in a
#! path of q points where the list of q points is gicen in 2pi/a units, with
#! a the lattice parameter given in Arnstrong

CC.Spectral.get_diag_dynamic_correction_along_path(dyn=dyn,
                                                   tensor3=tensor3,
                                                   k_grid=k_grid,
                                                   q_path_file="XGX.dat",
                                                   T = 200.0,
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   filename_sp = 'nomm_spectral_func_in_path')

## post processing:

1. w, p = minimizer.dyn.DiagonalizeSupercell()
2. relax.minim.plot_results()
   1. ens.update_weights(final_dyn, T)
3. free_energy_hessian = relax.minim.ensemble.get_free_energy_hessian(...)
