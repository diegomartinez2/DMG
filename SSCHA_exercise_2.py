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

class Calculo_inicial(object):
    def __init__(self,fichero_ForceFields,fichero_dyn,nqirr):
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
        self.minim.minim_struct = True

        self.relax = sscha.Relax.SSCHA(self.minim,
                          ase_calculator = self.ff_calculator,
                          N_configs = 1000,
                          max_pop = 200)

        # Setup the custom function to print the frequencies at each step of the minimization
        self.io_func = sscha.Utilities.IOInfo()
        self.io_func.SetupSaving(fichero_frecuencias) # The file that will contain the frequencies is frequencies.dat

        # Now tell relax to call the function to save the frequencies after each iteration
        # CFP stands for Custom Function Post (Post = after the minimization step)
        self.relax.setup_custom_functions(custom_function_post = self.io_func.CFP_SaveFrequencies)
        # Finalmente hacemos todos los calculos de busqueda de la energia libre.
        self.relax.relax()
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
    def __init__(self,fichero_ForceFields,fichero_dyn,nqirr):
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
        self.ensemble.get_energy_forces(self.ff_calculator, compute_stress = False) #test compute_stress = True no puede con este potencial...
    def hessiano(self,T):
        #self.dyn_hessian = self.ensemble.get_free_energy_hessian(include_v4 = False) # We neglect high-order four phonon scattering

        print("Updating the importance sampling...")
        self.ensemble.update_weights(self.dyn_sscha_final, T0 = T)

        print("Computing the free energy hessian...")
        self.dyn_hessian = self.ensemble.get_free_energy_hessian(include_v4 = INCLUDE_V4,
                                                  get_full_hessian = True,verbose = True) # Full calculus
        # We can save it
        self.dyn_hessian.save_qe("hessian")

        w_hessian, pols_hessian = self.dyn_hessian.DiagonalizeSupercell()

        # Print all the frequency converting them into cm-1 (They are in Ry)
        print("\n".join(["{:16.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in w_hessian]))

class Hessiano_Vs_Temperatura(object):
    def __init__(self,T0,temperatura_i,fichero_ForceFields):
        # Load the dynamical matrix for the force field
        self.ff_dyn = CC.Phonons.Phonons(fichero_ForceFields, 3)

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
            self.minim.min_step_struc = 0.05
            self.minim.min_step_dyn = 0.002
            self.minim.kong_liu_ratio = 0.5
            self.minim.meaningful_factor = 0.000001
            #minim.root_representation = "root4"
            #minim.precond_dyn = False
            self.minim.minim_struct = True # *test*
#            self.minim.neglect_symmetries = True    # *test*

            # Prepare the relaxer (through many population)
#            self.relax = sscha.Relax.SSCHA(self.minim, ase_calculator = self.ff_calculator, N_configs=1000, max_pop=50)
            self.relax = sscha.Relax.SSCHA(self.minim, ase_calculator = self.ff_calculator, N_configs=1000, max_pop=200)

            # Relax
            self.relax.relax()

            # Save the dynamical matrix
            self.relax.minim.dyn.save_qe(Fichero_final_matriz_dinamica.format(int(Temperatura)))

            # Detect space group
            symm=spglib.get_spacegroup(self.relax.minim.dyn.structure.get_ase_atoms(), 0.005)
            print('Current SG = ', symm,' at T=',int(Temperatura))

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
            #print ("\n".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in pols_hessian]))
            #exit()

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

class Funcion_espectral(object):
    def __init__(self,Fichero_dyn_SnTe,nqirr):
        self.dyn = CC.Phonons.Phonons(Fichero_dyn_SnTe,nqirr)
        self.supercell = self.dyn.GetSupercell()
    def prepara_tensor(self):
        self.tensor3 =  CC.ForceTensor.Tensor3(self.dyn.structure,
                                self.dyn.structure.generate_supercell(self.supercell),
                                self.supercell)
        #! Assign the tensor3 values
        d3 = np.load("d3_realspace_sym.npy")*2.0 # The 2 factor is because of units, needs to be passed to Ry
        self.tensor3.SetupFromTensor(d3)
        #! Center and apply ASR, which is needed to interpolate the third order force constant
        self.tensor3.Center()
        self.tensor3.Apply_ASR()

        #! Print the tensor if you want, uncommenting the next line
        #self.tensor3.WriteOnFile(fname="FC3",file_format='D3Q')
    def calcula_espectro1(self,T0):
        # integration grid
        k_grid=[4,4,4]

        # q points in 2pi/Angstrom
        list_of_q_points=[ [  0.0000000,  0.0000000,  0.0000000 ],
                           [ -0.0386763,  0.0386763, -0.0386763 ],
                           [  0.0773527, -0.0773527,  0.0773527 ],
                           [  0.0000000,  0.0773527,  0.0000000 ],
                           [  0.1160290, -0.0386763,  0.1160290 ],
                           [  0.0773527,  0.0000000,  0.0773527 ],
                           [  0.0000000, -0.1547054,  0.0000000 ],
                           [ -0.0773527, -0.1547054,  0.0000000 ]   ]


        CC.Spectral.get_static_correction_along_path(dyn=self.dyn,
                                             tensor3=self.tensor3,
                                             k_grid=k_grid,
                                             q_path=list_of_q_points,
                                             filename_st="v2_v2+d3static_freq.dat",
                                             T = T0,
                                             print_dyn = False) # set true to print the Hessian dynamical matrices
                                                                # for each q point
    def calcula_espectro2(self,T0):
        # integration grid
        k_grid=[20,20,20]


        CC.Spectral.get_static_correction_along_path(dyn=self.dyn,
                                             tensor3=self.tensor3,
                                             k_grid=k_grid,
                                             q_path_file="XGX_path.dat",
                                             filename_st="v2_v2+d3static_freq.dat",
                                             T = T0,
                                             print_dyn = False) # set true to print the Hessian dynamical matrices
                                                                # for each q point
    def calcula_espectro3(self,T0):
        # integration grid
        k_grid=[20,20,20]

        # X and G in 2pi/Angstrom
        points=[[-0.1525326,  0.0,  0.0],
                [0.0       ,  0.0,  0.0]      ]

        CC.Spectral.get_full_dynamic_correction_along_path(dyn=self.dyn,
                                                   tensor3=self.tensor3,
                                                   k_grid=k_grid,
                                                   e1=100, de=0.1, e0=0,     # energy grid
                                                   sm1=1.0, sm0=1.0,  nsm=1, # smearing values
                                                   sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                                   T = T0,
                                                   q_path=points,
                                                   static_limit = True, #static approximation
                                                   notransl = True,  # projects out the acoustic zone center modes
                                                   filename_sp='static_spectral_func')

    def calcula_espectro4(self,T0):
        # integration grid
        k_grid=[20,20,20]

        # q point
        G=[0.0,0.0,0.0]


        CC.Spectral.get_full_dynamic_correction_along_path(dyn=self.dyn,
                                           tensor3=self.tensor3,
                                           k_grid=k_grid,
                                           e1=145, de=0.1, e0=0,
                                           sm1=1, sm0=1,nsm=1,
                                           sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                           T = T0,
                                           q_path=G,
                                           notransl = True,
                                           filename_sp='full_spectral_func')

    def calcula_espectro5(self,T0):
        # integration grid
        k_grid=[20,20,20]

        #
        G=[0.0,0.0,0.0]

        CC.Spectral.get_diag_dynamic_correction_along_path(dyn=self.dyn,
                                                   tensor3=self.tensor3,
                                                   k_grid=k_grid,
                                                   q_path=G,
                                                   T = T0,
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                                   filename_sp = 'nomm_spectral_func')

    def calcula_espectro6(self,T0):
        # integration grid
        k_grid=[20,20,20]

        CC.Spectral.get_diag_dynamic_correction_along_path(dyn=self.dyn,
                                                   tensor3=self.tensor3,
                                                   k_grid=k_grid,
                                                   q_path_file="XGX_path.dat",
                                                   T = T0,
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                                   filename_sp = 'nomm_spectral_func')


    def calcula_espectro(self,T0):
        """
        #! Calculate the spectral function at Gamma in the no-mode mixing approximation
        #! keeping the energy dependence on the self-energy.
        """
        #! An interpolation grid needs to be used (and one needs to check convergence with
        #! respect to it considering different values of the smearing)
        #! interpolation grid
        k_grid=[20,20,20]

        #
        G=[0.0,0.0,0.0]

        CC.Spectral.get_diag_dynamic_correction_along_path(dyn=self.dyn,
                                                   tensor3=self.tensor3,
                                                   k_grid=k_grid,
                                                   q_path=G,
                                                   T = T0,                             # The temperature for the calculation
                                                   e1=145, de=0.1, e0=0,                # The energy grid in cm-1
                                                   sm1=1.0, nsm=1, sm0=1.0,             # The smearing \eta for the analytic continuation
                                                   sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                                   filename_sp = 'nomm_spectral_func')  # Output file name

        #! Now perform the calculation of the spectral function in a
        #! path of q points where the list of q points is gicen in 2pi/a units, with
        #! a the lattice parameter given in Arnstrong

        CC.Spectral.get_diag_dynamic_correction_along_path(dyn=self.dyn,
                                                   tensor3=self.tensor3,
                                                   k_grid=k_grid,
                                                   q_path_file="XGX.dat",
                                                   T = T0,
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                                   filename_sp = 'nomm_spectral_func_in_path')

def main(args):
    #La temperatura del primer calculo
    T0 = 0
    #Las temperaturas de los otros calculos
    Temperatura_i = np.linspace(100, 250, 10)
    #El fichero de la matrix dinámica para el campo de fuerzas (entrada)
    Fichero_ForceFields = "ffield_dynq"
    #El fichero de la matriz dinamica del sistema SnTe
    #Fichero_dyn_SnTe = "start_dyn_"
    #Fichero_dyn_SnTe = "ffield_dynq"
    Fichero_dyn_SnTe = "final_sscha_T100_"
    #y el numero de ficheros, relacionado con q mesh del quantum espresso (y a su vez relacionado con la supercelda)
    nqirr = 4
    #El fichero de las frecuencias (salida)
    Fichero_frecuencias = "frequencies.dat"
    #Los ficheros de la matriz dinamica (salida)
    Fichero_final_matriz_dinamica = "final_sscha_T{}_"
    #"final_sscha_T{}_".format(int(T))

    Calculo = Calculo_inicial(Fichero_ForceFields,Fichero_dyn_SnTe,nqirr)
    Calculo.ensambla(T0)
    #Calculo.minimiza(Fichero_frecuencias,Fichero_final_matriz_dinamica.format(int(T0)),nqirr)
    Calculo.minimiza(Fichero_frecuencias,Fichero_final_matriz_dinamica.format(int(T0)))
    Calculo.dibuja(Fichero_frecuencias)

    Inestable = Busca_inestabilidades(Fichero_ForceFields,Fichero_dyn_SnTe,nqirr)
    Inestable.load_dyn(Fichero_final_matriz_dinamica.format(int(T0)),nqirr)
    Inestable.ensambla(T0)
    Inestable.calcula1()
    Inestable.hessiano(T0)

    Espectro =  Funcion_espectral(Fichero_dyn_SnTe,nqirr)
    Espectro.prepara_tensor()
    Espectro.calcula_espectro(T0)

    #aqui se mete el bucle en temperaturas para crear la entrada de datos a Hessiano_Vs_Temperatura
    ##temperatura_i = np.linspace(50, 300, 6)
    #for Temperatura in Temperatura_i:
    #   Calculo.ensambla(Temperatura)
    #   Calculo.minimiza(Fichero_frecuencias,Fichero_final_matriz_dinamica.format(int(Temperatura)))

    HessianoVsTemperatura = Hessiano_Vs_Temperatura(T0,Temperatura_i,Fichero_ForceFields)
    HessianoVsTemperatura.ciclo_T(Fichero_final_matriz_dinamica,nqirr)
    HessianoVsTemperatura.dibuja()
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
