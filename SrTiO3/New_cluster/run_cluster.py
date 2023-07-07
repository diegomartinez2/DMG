import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import sscha.SchaMinimizer, sscha.Relax

# Import the two python scripts for the cluster and espresso configurations
import espresso_calculator
import cluster

# Generate an ensemble with 300 configurations
dyn = CC.Phonons.Phonons("start_sscha", 3)
ensemble = sscha.Ensemble.Ensemble(dyn, 300)

# Get the espresso and cluster configurations
espresso_config = espresso_calculator.get_calculator()
cluster_config = cluster.configure_cluster()

# Setup the minimizer
minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

# Setup the automatic relaxation
relax = sscha.Relax.SSCHA(minimizer, espresso_config,
        N_configs=10,
        max_pop=3,
        save_ensemble=True,
        cluster=cluster_config)

# Setup the IO to save the minimization data and the frequencies
ioinfo = sscha.Utilities.IOInfo()
ioinfo.SetupSaving("minimization_data")

# Activate the data saving in the minimization
relax.setup_custom_functions(custom_function_post=ioinfo.CFP_SaveAll)

# Perform the NVT simulation
relax.relax(get_stress=True)

# Save the data
relax.minim.finalize()
relax.minim.dyn.save_qe("final_dyn")
