from __future__ import print_function
from __future__ import division
# Import the modules to read the dynamical matrix
import cellconstructor as CC
import cellconstructor.Phonons
# Import the SCHA modules
import sscha, sscha.Ensemble

import numpy as np
import cellconstructor.ForceTensor
import cellconstructor.ThermalConductivity
import time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def thermal_calculo(d3, dyn_prefix = 'final_dyn',nqirr = 8,T_init=250,T_end=350,T_steps=3):
    """
    Calculate the lattice thermal conductivity with tho methods:
    * Boltzman transport equation (label as SRTA)
    * Green_Kubo method (label as GK)
    """
    SSCHA_TO_MS = cellconstructor.ThermalConductivity.SSCHA_TO_MS
    RY_TO_THZ = cellconstructor.ThermalConductivity.SSCHA_TO_THZ
    dyn = CC.Phonons.Phonons(dyn_prefix, nqirr)

    supercell = dyn.GetSupercell()

    fc3 = CC.ForceTensor.Tensor3(dyn.structure,
    dyn.structure.generate_supercell(supercell), supercell)

    #d3 = np.load("d3.npy")
    fc3.SetupFromTensor(d3)
    fc3 = CC.ThermalConductivity.centering_fc3(fc3)

    mesh = [10,10,10]
    smear = 0.03/RY_TO_THZ

    tc = CC.ThermalConductivity.ThermalConductivity(dyn, fc3,
    kpoint_grid = mesh, scattering_grid = mesh, smearing_scale = None,
    smearing_type = 'adaptive', cp_mode = 'quantum', off_diag = False) #'adaptive' or 'constant'

    temperatures = np.linspace(T_init,T_end,T_steps,dtype=float)
    start_time = time.time()
    tc.setup_harmonic_properties(smear)
    tc.write_harmonic_properties_to_file()

    tc.calculate_kappa(mode = 'SRTA', temperatures = temperatures,
    write_lifetimes = True, gauss_smearing = True, offdiag_mode = 'wigner',
    kappa_filename = 'Thermal_conductivity_SRTA', lf_method = 'fortran-LA')

    print('Calculated SSCHA kappa in: ', time.time() - start_time)

    tc.calculate_kappa(mode = 'GK', write_lineshapes=False,
    ne = 1000, temperatures = temperatures,     #ne=1000 ->ne=10000
    kappa_filename = 'Thermal_conductivity_GK')

    print('Calculated SSCHA kappa in: ', time.time() - start_time)
    # Save ThermalConductivity object for postprocessing.
    tc.save_pickle()

#----------loop for figure Tau vs. smear----------------------
def thermal_calculo_bis(d3, dyn_prefix = 'final_dyn',nqirr = 8,T=[300],Mesh_range=[10,30,10],smear_range=[0.01,0.1,10]):
    """
    Alternate function to calculate the lattice thermal conductance with two Methods
    within a ranged defined calculation mesh and smear values.
    """
    SSCHA_TO_MS = cellconstructor.ThermalConductivity.SSCHA_TO_MS
    RY_TO_THZ = cellconstructor.ThermalConductivity.SSCHA_TO_THZ
    dyn = CC.Phonons.Phonons(dyn_prefix, nqirr)

    supercell = dyn.GetSupercell()

    fc3 = CC.ForceTensor.Tensor3(dyn.structure,
    dyn.structure.generate_supercell(supercell), supercell)

    #d3 = np.load("d3.npy")
    fc3.SetupFromTensor(d3)
    fc3 = CC.ThermalConductivity.centering_fc3(fc3)

    for index_mesh in range(Mesh_range[0],Mesh_range[1],Mesh_range[2]):
        mesh = [index_mesh,index_mesh,index_mesh]
        for index_smear in np.linspace(smear_range[0],smear_range[1],smear_range[2]):
            smear = index_smear/RY_TO_THZ

            tc = CC.ThermalConductivity.ThermalConductivity(dyn, fc3,
            kpoint_grid = mesh, scattering_grid = mesh, smearing_scale = None,
            smearing_type = 'adaptive', cp_mode = 'quantum', off_diag = False) # 'adaptive' or 'constant'

            #temperatures = np.linspace(T_init,T_end,T_steps,dtype=float)
            start_time = time.time()
            tc.setup_harmonic_properties(smear)
            tc.write_harmonic_properties_to_file()

            tc.calculate_kappa(mode = 'SRTA', temperatures = T,
            write_lifetimes = True, gauss_smearing = True, offdiag_mode = 'wigner',
            kappa_filename = 'Thermal_conductivity_SRTA_Mesh{0}_smear{1}'.format(index_mesh,index_smear),
             lf_method = 'fortran-LA')

            print('Calculated SSCHA kappa in: ', time.time() - start_time)

            tc.calculate_kappa(mode = 'GK', write_lineshapes=False,
            ne = 2000, temperatures = T,     #ne=1000 ->ne=10000 Less anharmonic materials and lower temperatures will need more points (Number of frequency points to calculate phonon lineshapes in Gree-kubo).
            kappa_filename = 'Thermal_conductivity_GK_Mesh{0}_smear{1}'.format(index_mesh,index_smear))

            print('Calculated SSCHA kappa in: ', time.time() - start_time)
    # Save ThermalConductivity object for postprocessing.
    tc.save_pickle()
#----------------------------------------------------------------

def processing():
    """manual"""
    tc = CC.ThermalConductivity.load_thermal_conductivity()

    # See at which temperatures we calculated stuff
    tc.what_temperatures()

    key = list(tc.lineshapes.keys()) # Get Ts for lineshapes

    # DOS calculated from auxiliary force constants
    harm_dos = tc.get_dos()
    # Temperature dependent DOS calculated from lineshapes
    # first two arrays are raw data
    # second two is gaussian smoothed results \
    #for the distance between energy points de
    anharm_dos = tc.get_dos_from_lineshapes(float(key[-1]), de = 0.1)
    return harm_dos,anharm_dos,key

def plot(harm_dos, anharm_dos,key):
    """# Plot results"""
    fig = plt.figure(figsize=(6.4, 4.8))
    gs1 = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs1[0, 0])
    ax.plot(harm_dos[0], harm_dos[1], 'k-',
    zorder=0, label = 'Harmonic')
    ax.plot(anharm_dos[0], anharm_dos[1], 'r-',
    zorder=0, label = 'Anharmonic raw @ ' + key[-1] + ' K')
    ax.plot(anharm_dos[2], anharm_dos[3], 'b-',
    zorder=0, label = 'Anharmonic smooth @ ' + key[-1] + ' K')
    ax.set_xlabel('Frequency (THz)')
    ax.set_ylabel('Density of states')
    ax.legend(loc = 'upper right')
    ax.set_ylim(bottom = 0.0)
    fig.savefig('test.pdf')
    plt.show()
    pass

def Hessian_calculus(DATA_DIR = 'pop3/data',N_RANDOM = 512,DYN_PREFIX =  'pop3/dyn_start_population3_',
                    FINAL_DYN =   'pop3/dyn_end_population3_',SAVE_PREFIX = 'dyn_hessian_',
                    NQIRR = 4, Tg = 300, T =  300, POPULATION = 3, INCLUDE_V4 = False):
    # Here the input information
    #DATA_DIR = 'pop3/data'               # path to the directory ens_pop#lastpop where the last population is stored
    #N_RANDOM = 512                                  # number elements in the ensamble
    #DYN_PREFIX =  'pop3/dyn_start_population3_'          # dyn mat that generated the last population
    #FINAL_DYN =   'pop3/dyn_end_population3_'            # SSCHA dyn mat obtained with the last minimization
    #SAVE_PREFIX = 'dyn_hessian_'                    # Free energy Hessian dynamical matrices outpu
    #NQIRR = 4                                       # The number or irredubcible q points
    #Tg = 300                                         # The temperature used to generate the configurations
    #T =  300                                         # The temperature for the calculation
    #POPULATION = 3                                  # number of last population
    #INCLUDE_V4 = False                              # True to include the 4th-order SSCHA FC term to calculate the Hessian
    INFO = """
    In this example we compute the free energy hessian.

    The ensemble has been generated with the dynamical matrix at:
    {}

    And to compute the hessian we will use reweighting at:
    {}

    The original temperature was {} K, we use reweighting to {} K.
    The ensemble, population {}, is located at: {}
    The number of configuration is {}.
    Do we include the v4 in the calculation? {}

    The free energy Hessian will be saved in: {}

    The (symmetrized) 3rd order FCs in d3_realspace_sym.npy

    """.format(DYN_PREFIX, FINAL_DYN, Tg, T, POPULATION, DATA_DIR,
               N_RANDOM, INCLUDE_V4, SAVE_PREFIX)
    print(INFO)
    print()
    print(" ======== RUNNING ======== ")
    print()

    print("Loading the original dynamical matrix...")
    dyn = CC.Phonons.Phonons(DYN_PREFIX, NQIRR)
    print("Loading the current dynamical matrix...")
    final_dyn = CC.Phonons.Phonons(FINAL_DYN, NQIRR)

    print("Loading the ensemble...")
    ens = sscha.Ensemble.Ensemble(dyn, Tg, dyn.GetSupercell())
    ens.load(DATA_DIR, POPULATION, N_RANDOM)
    # If the ensemble was saved in binary format, load it with
    # ens.load_bin(DATA_DIR, POPULATION)

    print("Updating the importance sampling...")
    ens.update_weights(final_dyn, T)

    print("Computing the free energy hessian...")
    # Set get_full_hessian to false to have only the odd correction
    # Usefull if you want to study the convergence with the number of configuration
    dyn_hessian, d3  = ens.get_free_energy_hessian(include_v4 = INCLUDE_V4,
                                              get_full_hessian = True,
                                              verbose = False,
                                                  return_d3 = True)

    print("Saving the hessian to {}...".format(SAVE_PREFIX))
    dyn_hessian.save_qe(SAVE_PREFIX)
    print("Hessian done.")
    ############################### FC3 part #############################################

    tensor3 = CC.ForceTensor.Tensor3(dyn=dyn_hessian)          # initialize 3rd order tensor
    tensor3.SetupFromTensor(d3)                              # assign values
    tensor3.Center()                                         # center it
    tensor3.Apply_ASR()                                      # apply ASR
    tensor3.WriteOnFile(fname="FC3",file_format='D3Q')       # write on file
    print("Done.")
    return d3

def main(args):
    """
    Main code for the thermal calculation using SSCHA.
    This code takes arguments from the command line:
    Population: the polulation index for the hessian calculation
    nqirr: the number of irreducible q points, number of files in the dynamical matrix output
    Tg: temperature used to generate the SSCHA configurations
    T: temperature of the calculation
    n_random: number of configurations used in teh essemble used in SSCHA

    Note: SrTiO3 undergoes a structural transition at 105 K [H. Unoki and T. Sakudo, Journal of the Physical Society of Japan 23, 546 (1967), https://doi.org/10.1143/JPSJ.23.546.]

    """
    if (len( sys.argv ) > 1):

        N_RANDOM = arg[5]

        SAVE_PREFIX = 'dyn_hessian_'
        NQIRR = arg[2]
        Tg = arg[3]
        T =  arg[4]
        POPULATION = arg[1]
        DATA_DIR = "pop{}/data".format(POPULATION)
        DYN_PREFIX =  "pop{}/dyn_start_population3_".format(POPULATION)
        FINAL_DYN =   "pop{}/dyn_end_population3_".format(POPULATION)
        INCLUDE_V4 = False
        print ("Population:",POPULATION)
        print ("Number elements in the ensamble",N_RANDOM)
        print ("The number or irredubcible q points (nqirr):",NQIRR)
        print ("The temperature used to generate the configurations:",Tg)
        print ("The temperature for the calculation:",T)
        print ("Path to the directory ens_pop#lastpop where the last population is stored:",DATA_DIR)
        print ("The dynamical matrix that generated the last population:",DYN_PREFIX)
        print ("The SSCHA dynamical matrix obtained with the last minimization",FINAL_DYN)
        print ("Free energy Hessian dynamical matrices output:",SAVE_PREFIX)
        d3 = Hessian_calculus(DATA_DIR,N_RANDOM,DYN_PREFIX,FINAL_DYN,SAVE_PREFIX,
                        NQIRR,Tg,T,POPULATION,INCLUDE_V4)
        thermal_calculo(d3,FINAL_DYN,NQIRR)
        harm_dos, anharm_dos ,key = processing()
        plot(harm_dos, anharm_dos, key)
        np.savetxt("dos_harmonic.dat",harm_dos,header='Temperature dependent Harmonic DOS from auxiliary force constants:')
        np.savetxt("dos_anharmonic.dat",anharm_dos,header='Temperature dependent Anharmonic DOS from lineshapes: 2 lines of raw data 2 lines of gaussian smoothed data')
    else:
        print ("Arguments are population, nqirr, Tg, T ,nrandom. In that order, separated by simple spaces")
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
