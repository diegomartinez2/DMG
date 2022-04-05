#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SSCHA_exercise.py
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

# Import the cellconstructor stuff
import cellconstructor as CC
import cellconstructor.Phonons

# Import the modules of the force field
import fforces as ff
import fforces.Calculator

# Import the modules to run the sscha
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.Utilities

# Import Matplotlib to plot
import numpy as np
import matplotlib.pyplot as plt

class Calculo_inicial(object):
    def __init__(self,fichero_ForceFields,nqirr):
        # Load the dynamical matrix for the force field
        self.ff_dyn = CC.Phonons.Phonons(fichero_ForceFields, nqirr)

        # Setup the forcefield with the correct parameters
        self.ff_calculator = ff.Calculator.ToyModelCalculator(self.ff_dyn)
        self.ff_calculator.type_cal = "pbtex"
        self.ff_calculator.p3 = 0.036475
        self.ff_calculator.p4 = -0.022
        self.ff_calculator.p4x = -0.014

        # Initialization of the SSCHA matrix
        self.dyn_sscha = self.ff_dyn.Copy()
        # Flip the imaginary frequencies into real ones
        self.dyn_sscha.ForcePositiveDefinite()
        # Apply the ASR and the symmetry group
        self.dyn_sscha.Symmetrize()

    def ensambla(self,T):
        self.ensemble = sscha.Ensemble.Ensemble(self.dyn_sscha, T0 = T, supercell = self.dyn_sscha.GetSupercell())

    def minimiza(self,fichero_frecuencias,fichero_matriz):
        self.minim = sscha.SchaMinimizer.SSCHA_Minimizer(self.ensemble)

        # Lets setup the minimization on the fourth root
        self.minim.root_representation = "root4" # Other possibilities are 'normal' and 'sqrt'

        # To work correctly with the root4, we must deactivate the preconditioning on the dynamical matrix
        self.minim.precond_dyn = False

        # Now we setup the minimization parameters
        # Since we are quite far from the correct solution, we will use a small optimization step
        self.minim.min_step_dyn = 1 # If the minimization ends with few steps (less than 10), decrease it, if it takes too much, increase it

        # We decrease the Kong-Liu effective sample size below which the population is stopped
        self.minim.kong_liu_ratio = 0.2 # Default 0.5

        self.relax = sscha.Relax.SSCHA(self.minim,
                          ase_calculator = self.ff_calculator,
                          N_configs = 1000,
                          max_pop = 20)

        # Setup the custom function to print the frequencies at each step of the minimization
        self.io_func = sscha.Utilities.IOInfo()
        self.io_func.SetupSaving(fichero_frecuencias) # The file that will contain the frequencies is frequencies.dat

        # Now tell relax to call the function to save the frequencies after each iteration
        # CFP stands for Custom Function Post (Post = after the minimization step)
        self.relax.setup_custom_functions(custom_function_post = self.io_func.CFP_SaveFrequencies)
        # Finalmente hacemos todos los calculos de busqueda de la energia libre.
        self.relax.relax()

        # Save the final dynamical matrix
        self.relax.minim.dyn.save_qe(fichero_matriz)

    def dibuja(self,fichero):
        # Setup the interactive plotting mode
        #plt.ion()

        # Lets plot the Free energy, gradient and the Kong-Liu effective sample size
        self.relax.minim.plot_results()

        frequencies = np.loadtxt(fichero)
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
        plt.savefig('Step_Freq.png')

class Busca_inestabilidades(object):
    def __init__(self,fichero_ForceFields,nqirr):
        # Load the dynamical matrix for the force field
        self.ff_dyn = CC.Phonons.Phonons(fichero_ForceFields, nqirr)

        # Setup the forcefield with the correct parameters
        self.ff_calculator = ff.Calculator.ToyModelCalculator(self.ff_dyn)
        self.ff_calculator.type_cal = "pbtex"
        self.ff_calculator.p3 = 0.036475
        self.ff_calculator.p4 = -0.022
        self.ff_calculator.p4x = -0.014

        # Initialization of the SSCHA matrix
        self.dyn_sscha = self.ff_dyn.Copy()
        self.dyn_sscha.ForcePositiveDefinite()

        # Apply also the ASR and the symmetry group
        self.dyn_sscha.Symmetrize()
    def load_dyn(self,Fichero_final_matriz_dinamica,nqirr):
        # We reload the final result (no need to rerun the sscha minimization)
        self.dyn_sscha_final = CC.Phonons.Phonons(Fichero_final_matriz_dinamica, nqirr)
    def ensambla(self,T):
        # We reset the ensemble
        self.ensemble = sscha.Ensemble.Ensemble(self.dyn_sscha_final, T0 = T, supercell = self.dyn_sscha_final.GetSupercell())

        # We need a bigger ensemble to properly compute the hessian
        # Here we will use 10000 configurations
        self.ensemble.generate(10000)
    def calcula1(self):
        # We now compute forces and energies using the force field calculator
        self.ensemble.get_energy_forces(self.ff_calculator, compute_stress = False)
    def hessiano(self):
        self.dyn_hessian = self.ensemble.get_free_energy_hessian(include_v4 = False) # We neglect high-order four phonon scattering
        # We can save it
        self.dyn_hessian.save_qe("hessian")

        w_hessian, pols_hessian = self.dyn_hessian.DiagonalizeSupercell()

        # Print all the frequency converting them into cm-1 (They are in Ry)
        print("\n".join(["{:16.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in w_hessian]))

class Hessiano_Vs_Temperatura(object):
    def __init__(self,T0,temperatura_i,fichero_ForceFields,nqirr):
        # Load the dynamical matrix for the force field
        self.ff_dyn = CC.Phonons.Phonons(fichero_ForceFields, nqirr)

        # Setup the forcefield with the correct parameters
        self.ff_calculator = ff.Calculator.ToyModelCalculator(self.ff_dyn)
        self.ff_calculator.type_cal = "pbtex"
        self.ff_calculator.p3 = 0.036475
        self.ff_calculator.p4 = -0.022
        self.ff_calculator.p4x = -0.014
        # Define the temperatures, from 50 to 300 K, 6 temperatures
        #self.temperatures = np.linspace(50, 300, 6)
        self.temperatures = temperatura_i

        self.lowest_hessian_mode = []
        self.lowest_sscha_mode = []

        # Perform a simulation at each temperature
        self.t_old = T0

    def ciclo_T(self,Fichero_final_matriz_dinamica,nqirr):
        for Temperatura in self.temperatures:
            # Load the starting dynamical matrix
            self.dyn = CC.Phonons.Phonons(Fichero_final_matriz_dinamica.format(int(self.t_old)), nqirr)

            # Prepare the ensemble
            self.ensemble = sscha.Ensemble.Ensemble(self.dyn, T0 = Temperatura, supercell = self.dyn.GetSupercell())

            # Prepare the minimizer
            self.minim = sscha.SchaMinimizer.SSCHA_Minimizer(self.ensemble)
            self.minim.min_step_dyn = 0.002
            self.minim.kong_liu_ratio = 0.5
            #minim.root_representation = "root4"
            #minim.precond_dyn = False

            # Prepare the relaxer (through many population)
            self.relax = sscha.Relax.SSCHA(self.minim, ase_calculator = self.ff_calculator, N_configs=1000, max_pop=5)

            # Relax
            self.relax.relax()

            # Save the dynamical matrix
            self.relax.minim.dyn.save_qe(Fichero_final_matriz_dinamica.format(int(Temperatura)))

            # Recompute the ensemble for the hessian calculation
            self.ensemble = sscha.Ensemble.Ensemble(self.relax.minim.dyn, T0 = Temperatura, supercell = self.dyn.GetSupercell())
            self.ensemble.generate(5000)
            self.ensemble.get_energy_forces(self.ff_calculator, compute_stress = False) #gets the energies and forces from ff_calculator

            #update weights!!! es posible que este sea el motivo por el que no obtengo buenos resultados?
            self.ensemble.update_weights(self.relax.minim.dyn, Temperatura)
            # Get the free energy hessian
            dyn_hessian = self.ensemble.get_free_energy_hessian(include_v4 = False) #free energy hessian as in Bianco paper 2017
            dyn_hessian.save_qe("hessian_T{}_".format(int(Temperatura)))

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

            self.t_old = Temperatura
        # We prepare now the file to save the results
        freq_data = np.zeros( (len(self.temperatures), 3))
        freq_data[:, 0] = self.temperatures
        freq_data[:, 1] = self.lowest_sscha_mode
        freq_data[:, 2] = self.lowest_hessian_mode

        # Save results on file
        np.savetxt("hessian_vs_temperature.dat", freq_data, header = "T [K]; SSCHA mode [cm-1]; Free energy hessian [cm-1]")

    def dibuja(self):
        hessian_data = np.loadtxt("hessian_vs_temperature.dat")

        plt.figure(dpi = 120)
        plt.plot(hessian_data[:,0], hessian_data[:,1], label = "Min SCHA freq", marker = ">")
        plt.plot(hessian_data[:,0], hessian_data[:,2], label = "Free energy curvature", marker = "o")
        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("Temperature [K]")
        plt.ylabel("Frequency [cm-1]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('Temp_Freq.png')
        #plt.show()

        plt.figure(dpi = 120)
        plt.plot(hessian_data[:,0], np.sign(hessian_data[:,2]) * hessian_data[:,2]**2, label = "Free energy curvature", marker = "o")
        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("Temperature [K]")
        plt.ylabel("$\omega^2$ [cm-2]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('Temp_Omeg.png')
        #plt.show()
class Calculo_general(object):
    def __init__(self,T0,temperatura_i,fichero_ForceFields,nqirr):
        #dyn = CC.Phonons.Phonons(namefile, NQIRR)
        # Load the dynamical matrix for the force field
        self.ff_dyn = CC.Phonons.Phonons(fichero_ForceFields, nqirr)

        # Setup the forcefield with the correct parameters
        self.ff_calculator = ff.Calculator.ToyModelCalculator(self.ff_dyn)
        self.ff_calculator.type_cal = "pbtex"
        self.ff_calculator.p3 = 0.036475
        self.ff_calculator.p4 = -0.022
        self.ff_calculator.p4x = -0.014
        self.dyn = self.ff_dyn.Copy()

        self.dyn.Symmetrize()

        self.dyn.ForcePositiveDefinite()
    def ensambla(self,Temperatura):
        self.ensemble = sscha.Ensemble.Ensemble(self.dyn, T0 = Temperatura, supercell = self.dyn.GetSupercell())
        #ensemble = sscha.Ensemble.Ensemble(dyn, 1000, supercell = dyn.GetSupercell())
        #ensemble.generate(N)
    def minimiza(self):
        self.minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(self.ensemble)
        self.minimizer.min_step_dyn = 0.005         # The minimization step on the dynamical matrix
        self.minimizer.min_step_struc = 0.05        # The minimization step on the structure
        self.minimizer.kong_liu_ratio = 0.5         # The parameter that estimates whether the ensemble is still good
        self.minimizer.meaningful_factor = 0.000001 # How much small the gradient should be before I stop?
    def paso_a_paso(self):
        self.minimizer.init()
        self.minimizer.run()
        #...
        self.minimizer.finalize()
        self.minimizer.plot_results()
    def relaja(self):
        relax = sscha.Relax.SSCHA(self.minimizer,
                          ase_calculator = self.ff_calculator,
                          N_configs = 10000,
                          max_pop = 20)
        relax.relax()
    def hessiano(self,DYN_PREFIX,file_FINAL_DYN,NQIRR,Tg,DATA_DIR, POPULATION, N_RANDOM,SAVE_PREFIX):
        #* Matriz dinamica original:
        dyn = CC.Phonons.Phonons(DYN_PREFIX, NQIRR)

        #* Matriz dinamica actual:
        final_dyn = CC.Phonons.Phonons(file_FINAL_DYN, NQIRR)

        #* reensmblado:
        ens = sscha.Ensemble.Ensemble(dyn, Tg, dyn.GetSupercell())
        ens.load(DATA_DIR, POPULATION, N_RANDOM)

        #* Importante reajustar los pesos:
        ens.update_weights(final_dyn, Tg)

        dyn_hessian = ens.get_free_energy_hessian(include_v4 = INCLUDE_V4,
                                          get_full_hessian = True,
                                          verbose = True)
        dyn_hessian.save_qe(SAVE_PREFIX)
    def spectral_function(self):
        # Initialize the tensor3 object
        # We need 2nd FCs of the used grid to configure the supercell.
        # For example, we can use the sscha final auxiliary dynamical matrices
        dyn = CC.Phonons.Phonons("dyn_end_population3_",3)
        supercell = dyn.GetSupercell()
        tensor3 = CC.ForceTensor.Tensor3(dyn.structure,
                                dyn.structure.generate_supercell(supercell),
                                supercell)

       # Assign the tensor3 values
       d3 = np.load("d3_realspace_sym.npy")*2.0 # The 2 factor is because of units, needs to be passed to Ry
       tensor3.SetupFromTensor(d3)

       # Center and apply ASR, which is needed to interpolate the third order force constant
       tensor3.Center()
       tensor3.Apply_ASR()

       # Print the tensor if you want, uncommenting the next line
       #tensor3.WriteOnFile(fname="FC3",file_format='D3Q')

       # Calculate the spectral function at Gamma in the no-mode mixing approximation
       # keeping the energy dependence on the self-energy.
       #
       # An interpolation grid needs to be used (and one needs to check convergence with
       # respect to it considering different values of the smearing)


       # interpolation grid
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

       # Now perform the calculation of the spectral function in a
       # path of q points where the list of q points is gicen in 2pi/a units, with
       # a the lattice parameter given in Arnstrong

       CC.Spectral.get_diag_dynamic_correction_along_path(dyn=dyn,
                                                   tensor3=tensor3,
                                                   k_grid=k_grid,
                                                   q_path_file="XGX.dat",
                                                   T = 200.0,
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   filename_sp = 'nomm_spectral_func_in_path')


def main(args):
    #La temperatura del primer calculo
    T0 = 0
    #Las temperaturas de los otros calculos
    Temperatura_i = np.linspace(50, 300, 6)
    #El fichero de la matrix dinámica para el campo de fuerzas (entrada)
    Fichero_ForceFields = "ffield_dynq"
    #y el numero de ficheros, relacionado con q mesh del quantum espresso (y a su vez relacionado con la supercelda)
    nqirr = 3
    #El fichero de las frecuencias (salida)
    Fichero_frecuencias = "frequencies.dat"
    #Los ficheros de la matriz dinamica (salida)
    Fichero_final_matriz_dinamica = "final_sscha_T{}_"
    #"final_sscha_T{}_".format(int(T))

    Calculo = Calculo_inicial(Fichero_ForceFields,nqirr)
    Calculo.ensambla(T0)
    Calculo.minimiza(Fichero_frecuencias,Fichero_final_matriz_dinamica.format(int(T0)),nqirr)
    Calculo.dibuja(Fichero_frecuencias)

    Inestable = Busca_inestabilidades(Fichero_ForceFields,nqirr)
    Inestable.load_dyn(Fichero_final_matriz_dinamica.format(int(T0),nqirr))
    Inestable.ensambla(T0)
    Inestable.calcula1()
    Inestable.hessiano()

    #aqui se mete el bucle en temperaturas para crear la entrada de datos a Hessiano_Vs_Temperatura
    ##temperatura_i = np.linspace(50, 300, 6)
    #for Temperatura in Temperatura_i:
    #   Calculo.ensambla(Temperatura)
    #   Calculo.minimiza(Fichero_frecuencias,Fichero_final_matriz_dinamica.format(int(Temperatura)))

    HessianoVsTemperatura = Hessiano_Vs_Temperatura(T0,Temperatura_i,Fichero_ForceFields,nqirr)
    HessianoVsTemperatura.ciclo_T(Fichero_final_matriz_dinamica,nqirr)
    HessianoVsTemperatura.dibuja()
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
