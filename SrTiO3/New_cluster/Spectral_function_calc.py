#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Spectral_function_calc.py
#
#  Copyright 2023 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
# ---------------------------
# Importación de los módulos
# ---------------------------
# Import the cellconstructor stuff
import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

# Import the modules of the force field
#import fforces as ff
#import fforces.Calculator

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
# -------
# Clases
# -------
class Funcion_espectral(object):
    """
    Object:
    """
    def __init__(self,Fichero_dyn,nqirr, PATH = "GXMGRX"):
        self.dyn = CC.Phonons.Phonons(Fichero_dyn,nqirr)
        self.supercell = self.dyn.GetSupercell()
        N_POINTS = 1000
        SPECIAL_POINTS = {"G": [0,0,0],
                      "X": [0, 0, .5],
                      "M": [0, .5, .5],
                      "R": [.5, .5, .5]}
        self.qpath, self.data = CC.Methods.get_bandpath(self.dyn.structure.unit_cell,
                                          PATH,
                                          SPECIAL_POINTS,
                                          N_POINTS)
    def prepara_tensor(self):
        self.tensor3 =  CC.ForceTensor.Tensor3(self.dyn.structure,
                                self.dyn.structure.generate_supercell(self.supercell),
                                self.supercell)
         #! Assign the tensor3 values
        d3 = np.load("d3_realspace_sym.npy")
        self.tensor3.SetupFromTensor(d3)
          #! Center and apply ASR, which is needed to interpolate the third order force constant
#        self.tensor3.Center()
#        self.tensor3.Apply_ASR()
        self.tensor3.Center(Far=2)
        self.tensor3.Apply_ASR(PBC=True)

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
    def dibuja1(self):
        plot_data = np.loadtxt("v2_v2+d3static_freq.dat")

        plt.figure(dpi = 120)
        plt.plot(plot_data[:,0], plot_data[:,1], marker = "o")
        plt.plot(plot_data[:,0], plot_data[:,2], marker = "o")
        plt.plot(plot_data[:,0], plot_data[:,3], marker = "o")
        plt.plot(plot_data[:,0], plot_data[:,4], marker = "o")
        plt.plot(plot_data[:,0], plot_data[:,5], marker = "o")
        plt.plot(plot_data[:,0], plot_data[:,6], marker = "o")
        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("Path (2pi/Angstrom)")
        plt.ylabel("Frequency [cm-1]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('v2_v2_path_Freq.png')
        #plt.show()

    def calcula_espectro2(self,T0):
        # integration grid
        k_grid=[20,20,20]


        CC.Spectral.get_static_correction_along_path(dyn=self.dyn,
                                             tensor3=self.tensor3,
                                             k_grid=k_grid,
                                             q_path_file="XGX_path.dat",
                                             filename_st="v2_v2+d3static_freq2.dat",
                                             T = T0,
                                             print_dyn = False) # set true to print the Hessian dynamical matrices
                                                                # for each q point
    def calcula_espectro2multiprocessing(self,T0,processes):
        # integration grid
        k_grid=[20,20,20]


        CC.Spectral.get_static_correction_along_path_multiprocessing(dyn=self.dyn,
                                             tensor3=self.tensor3,
                                             k_grid=k_grid,
                                             q_path_file="XGX_path.dat",
                                             filename_st="v2_v2+d3static_freq2_multiprocessing.dat",
                                             T = T0,
                                             print_dyn = False, processes = processes) # set true to print the Hessian dynamical matrices
                                                                # for each q point
    def dibuja2(self):
        plot_data = np.loadtxt("v2_v2+d3static_freq2.dat")

        plt.figure(dpi = 120)
        plt.plot(plot_data[:,0], plot_data[:,1])
        plt.plot(plot_data[:,0], plot_data[:,2])
        plt.plot(plot_data[:,0], plot_data[:,3])
        plt.plot(plot_data[:,0], plot_data[:,4])
        plt.plot(plot_data[:,0], plot_data[:,5])
        plt.plot(plot_data[:,0], plot_data[:,6])
#        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("XGX")
        plt.ylabel("Frequency [cm-1]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('v2_v2+d3static_freq2.png')
        #plt.show()
    def dibuja2multiprocessing(self):
        plot_data = np.loadtxt("v2_v2+d3static_freq2_multiprocessing.dat")

        plt.figure(dpi = 120)
        plt.plot(plot_data[:,0], plot_data[:,1])
        plt.plot(plot_data[:,0], plot_data[:,2])
        plt.plot(plot_data[:,0], plot_data[:,3])
        plt.plot(plot_data[:,0], plot_data[:,4])
        plt.plot(plot_data[:,0], plot_data[:,5])
        plt.plot(plot_data[:,0], plot_data[:,6])
#        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("XGX")
        plt.ylabel("Frequency [cm-1]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('v2_v2+d3static_freq2_multiprocessing.png')
        #plt.show()

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

    def dibuja3(self):
        plot_data = np.loadtxt("static_spectral_func_1.00_1.0.dat")

        plt.figure(dpi = 120)
        plt.plot(plot_data[:,1], plot_data[:,2])
#        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("Energy [cm-1]")
        plt.ylabel("Spectral Function [cm-1]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('static_spectral_func_1.00_1.0.png')
        #plt.show()

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
    def dibuja4(self):
        plot_data = np.loadtxt("full_spectral_func_1.00_1.0.dat")

        plt.figure(dpi = 120)
        plt.plot(plot_data[:,1], plot_data[:,2])
#        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("Energy [cm-1]")
        plt.ylabel("Spectral Function [cm-1]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('full_spectral_func_1.00_1.0.png')
        #plt.show()

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
    def dibuja5(self):
        plot_data = np.loadtxt("nomm_spectral_func_1.00.dat")

        plt.figure(dpi = 120)
        plt.plot(plot_data[:,1], plot_data[:,2])
        plt.plot(plot_data[:,1], plot_data[:,3])
        plt.plot(plot_data[:,1], plot_data[:,4])
        plt.plot(plot_data[:,1], plot_data[:,5])
        plt.plot(plot_data[:,1], plot_data[:,6])
        plt.plot(plot_data[:,1], plot_data[:,7])
        plt.plot(plot_data[:,1], plot_data[:,8])
#        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("Energy [cm-1]")
        plt.ylabel("Spectral Function [1/cm-1]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('nomm_spectral_func_1.00.png')
        #plt.show()

        plot_data = np.loadtxt("nomm_spectral_func_lorentz_one_shot_1.00.dat")

        plt.figure(dpi = 120)
        plt.plot(plot_data[:,1], plot_data[:,2])
        plt.plot(plot_data[:,1], plot_data[:,3])
        plt.plot(plot_data[:,1], plot_data[:,4])
        plt.plot(plot_data[:,1], plot_data[:,5])
        plt.plot(plot_data[:,1], plot_data[:,6])
        plt.plot(plot_data[:,1], plot_data[:,7])
        plt.plot(plot_data[:,1], plot_data[:,8])
#        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("Energy [cm-1]")
        plt.ylabel("Spectral Function [1/cm-1]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('nomm_spectral_func_lorentz_one_shot_1.00.png')
        #plt.show()

        plot_data = np.loadtxt("nomm_spectral_func_lorentz_perturb_1.00.dat")

        plt.figure(dpi = 120)
        plt.plot(plot_data[:,1], plot_data[:,2])
        plt.plot(plot_data[:,1], plot_data[:,3])
        plt.plot(plot_data[:,1], plot_data[:,4])
        plt.plot(plot_data[:,1], plot_data[:,5])
        plt.plot(plot_data[:,1], plot_data[:,6])
        plt.plot(plot_data[:,1], plot_data[:,7])
        plt.plot(plot_data[:,1], plot_data[:,8])
#        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("Energy [cm-1]")
        plt.ylabel("Spectral Function [1/cm-1]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('nomm_spectral_func_lorentz_perturb_1.00.png')
        #plt.show()


    def calcula_espectro6(self,T0):
        # integration grid
        k_grid=[20,20,20]

        CC.Spectral.get_diag_dynamic_correction_along_path(dyn=self.dyn,
                                                   tensor3=self.tensor3,
                                                   k_grid=k_grid,
                                                   q_path_file="XGX_path.dat",
                                                   T =T0,
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                                   filename_sp = 'nomm_spectral_func2')
    def dibuja6(self):
        # Prepare plot of phonon spectra
        fig, ax1 = plt.subplots(1,1)
        data = np.loadtxt('nomm_spectral_func2_1.00.dat')
        plt.scatter(data[:,0], data[:,1], s=1, c=data[:,2], cmap='hot')
        ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)
        plt.colorbar()
        plt.savefig('spectral_path.pdf', bbox_inches='tight')
        plt.savefig('nomm_spectral_func2_1.00.png')
        return 0

    def calcula_espectro6multiprocessing(self,T0,processes):
        # integration grid
        k_grid=[20,20,20]

        CC.Spectral.get_diag_dynamic_correction_along_path_multiprocessing(dyn=self.dyn,
                                                   tensor3=self.tensor3,
                                                   k_grid=k_grid,
                                                   q_path_file="XGX_path.dat",
                                                   T =T0,
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                                   filename_sp = 'nomm_spectral_func2_multiprocessing', processes = processes)
    def calcula_espectro6multiprocessing2(self,T0,processes):
        # integration grid
        k_grid=[20,20,20]

        CC.Spectral.get_diag_dynamic_correction_along_path_multiprocessing2(dyn=self.dyn,
                                                   tensor3=self.tensor3,
                                                   k_grid=k_grid,
                                                   q_path_file="XGX_path3.dat",
                                                   T =T0,
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                                   filename_sp = 'nomm_spectral_func2_multiprocessing2', processes = processes)
    def calcula_espectro6multiprocessing3(self,T0,processes):
        # integration grid
        k_grid=[20,20,20]

        CC.Spectral.get_diag_dynamic_correction_along_path_multiprocessing3(dyn=self.dyn,
                                                   tensor3=self.tensor3,
                                                   k_grid=k_grid,
                                                   q_path_file="XGX_path3.dat",
                                                   T =T0,
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                                   filename_sp = 'nomm_spectral_func2_multiprocessing3', processes = processes)
    def dibuja6multiprocessing(self):
        # Prepare plot of phonon spectra
        fig, ax1 = plt.subplots(1,1)
        data = np.loadtxt('nomm_spectral_func2_multiprocessing_1.00.dat')
        plt.scatter(data[:,0], data[:,1], s=1, c=data[:,2], cmap='hot')
        ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)
        plt.colorbar()
        plt.savefig('spectral_path2.pdf', bbox_inches='tight')
        plt.savefig('nomm_spectral_func2_multiprocessing_1.00.png')
        return 0
    def dibuja6multiprocessing2(self):
        # Prepare plot of phonon spectra
        fig, ax1 = plt.subplots(1,1)
        data = np.loadtxt('nomm_spectral_func2_multiprocessing_1.00.dat')
        X = data[:,0]
        Y = data[:,1]
        Z = data[:,2]
        z = [Z[i] for i in np.lexsort((Y,X))]
        data1 = np.resize(z,(int(len(X)/1450),1450))
        cax = ax1.imshow(data1.transpose(), cmap=cm.coolwarm, origin='lower')
        ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)
        xaxis, xticks, xlabels = self.data # Info to plot correclty the x axis
        plt.xticks(ticks=[0,499,999], labels=['X','G','X'])
        plt.yticks(ticks=[0,200,400,600,800,1000], labels=['0','20','40','60','80','100'])
        cbar = fig.colorbar(cax)
        plt.savefig('spectral_path2_imshow.pdf', bbox_inches='tight')
        plt.savefig('nomm_spectral_func2_multiprocessing_imshow_1.00.png')
        plt.show()
        return 0

    def calcula_espectro_full_Gamma(self,T0):
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
                                           filename_sp='full_spectral_Gamma')

    def calcula_espectro_full_R(self,T0):
        # integration grid
        k_grid=[20,20,20]

        # q point
        G=[0.5,0.5,0.5]


        CC.Spectral.get_full_dynamic_correction_along_path(dyn=self.dyn,
                                           tensor3=self.tensor3,
                                           k_grid=k_grid,
                                           e1=145, de=0.1, e0=0,
                                           sm1=1, sm0=1,nsm=1,
                                           sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                           T = T0,
                                           q_path=G,
                                           notransl = True,
                                           filename_sp='full_spectral_R')


#-------------------------SrTiO3--------------
    def calcula_espectro_basico_SrTiO3(self,T0):
        # integration grid
        k_grid=[20,20,20]


        CC.Spectral.get_static_correction_along_path(dyn=self.dyn,
                                             tensor3=self.tensor3,
                                             k_grid=k_grid,
                                             q_path=self.qpath,
                                             filename_st="SrTiO3_static.dat",
                                             T = T0,
                                             print_dyn = False) # set true to print the Hessian dynamical matrices
                                                                # for each q point
    def calcula_espectro_basico_SrTiO3_multiprocessing(self,T0, processes = None):
        # integration grid
        k_grid=[20,20,20]


        CC.Spectral.get_static_correction_along_path_multiprocessing(dyn=self.dyn,
                                             tensor3=self.tensor3,
                                             k_grid=k_grid,
                                             q_path=self.qpath,
                                             filename_st="SrTiO3_static_multiprocessing.dat",
                                             T = T0,
                                             print_dyn = False, processes = processes) # set true to print the Hessian dynamical matrices
                                                                # for each q point
    def dibuja_espectro_basico_SrTiO3(self,filename = "SrTiO3_static.dat", PATH = "GXMGRX"):
        plot_data = np.loadtxt(filename)

        nmodes = self.dyn.structure.N_atoms * 3
        # Plot the two dispersions
        plt.figure(dpi = 150)
        ax = plt.gca()
        xaxis, xticks, xlabels = self.data # Info to plot correclty the x axis
        for i in range(nmodes):
            lbl=None
            lblsscha = None
            if i == 0:
                lbl = 'Static'
                lblsscha = 'Static+bubble'

            ax.plot(xaxis, plot_data[:,i+1]*0.124, color = 'k', ls = 'dashed', label = lbl)
            ax.plot(xaxis, plot_data[:,i+nmodes+1]*0.124, color = 'r', label = lblsscha)

        # Plot vertical lines for each high symmetry points
        for x in xticks:
            ax.axvline(x, 0, 1, color = "k", lw = 0.4)
        ax.axhline(0, 0, 1, color = 'k', ls = ':', lw = 0.4)

        ax.legend()

        # Set the x labels to the high symmetry points
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels)

        ax.set_xlabel("Q path")
    #    ax.set_ylabel("Phonons [cm-1]")
        ax.set_ylabel("Phonons [meV]")

        plt.tight_layout()
        plt.savefig("{}_dispersion{}.png".format(filename, PATH))
        #plt.show()

    def calcula_espectro_correction_SrTiO3(self,T0):
        # integration grid
        k_grid=[20,20,20]

        CC.Spectral.get_diag_dynamic_correction_along_path(dyn=self.dyn,
                                                   tensor3=self.tensor3,
                                                   k_grid=k_grid,
                                                   q_path=self.qpath,
                                                   T =T0,
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                                   filename_sp = 'nomm_spectral_func2')
    def calcula_espectro_correction_multiprocessing_SrTiO3(self,T0, processes = None):
        # integration grid
        k_grid=[20,20,20]

        CC.Spectral.get_diag_dynamic_correction_along_path_multiprocessing(dyn=self.dyn,
                                                   tensor3=self.tensor3,
                                                   k_grid=k_grid,
                                                   q_path=self.qpath,
                                                   T =T0,
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   sm1_id=1.0, sm0_id=1.0,   # Minimum and maximum value of the smearing (cm-1) for the term of the Green function proportional to the identity
                                                   filename_sp = 'nomm_spectral_func2', processes = processes)

    def calcula_full_correction_en_punto_G(self,T0, processes = None):
        # integration grid
        k_grid=[20,20,20]

        CC.Spectral.get_full_dynamic_correction_along_path_multiprocessing(dyn=self.dyn,
                                           tensor3=self.tensor3,
                                           k_grid=k_grid,
                                           e1=100, de=0.1, e0=0,     # energy grid
                                           sm1=10.0,sm0=1.0,nsm=3,    # smearing values
                                           T=300,
                                           q_path=[0.0,0.0,0.0],
                                           static_limit = True,
                                           filename_sp='full_spectral_func_X', processes = processes)

    def calcula_full_correction_en_punto_R(self,T0, processes = None):
        # integration grid
        k_grid=[20,20,20]

        CC.Spectral.get_full_dynamic_correction_along_path_multiprocessing(dyn=self.dyn,
                                           tensor3=self.tensor3,
                                           k_grid=k_grid,
                                           e1=100, de=0.1, e0=0,     # energy grid
                                           sm1=10.0,sm0=1.0,nsm=3,    # smearing values
                                           T=300,
                                           q_path=[0.5,0.5,0.5],
                                           static_limit = True,
                                           filename_sp='full_spectral_func_X', processes = processes)


    def calcula_oneshot_correction_en_punto_Gamma(self,T0):


        # integration grid
        k_grid=[20,20,20]

        CC.Spectral.get_os_perturb_dynamic_correction_along_path(self.dyn,
                                           tensor3=self.tensor3,
                                           k_grid=k_grid,
                                           sm1=1.0, sm0=1.0,
                                           nsm=1,
                                           q_path=[0.0,0.0,0.0],
                                           T=T0)
    def calcula_oneshot_correction_en_punto_R(self,T0):


        # integration grid
        k_grid=[20,20,20]

        CC.Spectral.get_os_perturb_dynamic_correction_along_path(self.dyn,
                                           tensor3=self.tensor3,
                                           k_grid=k_grid,
                                           sm1=1.0, sm0=1.0,
                                           nsm=1,
                                           q_path=[0.5,0.5,0.5],
                                           T=T0)

    def calcula_oneshot_correction_along_PATH(self,T0):


        # integration grid
        k_grid=[20,20,20]

        CC.Spectral.get_os_perturb_dynamic_correction_along_path(self.dyn,
                                           tensor3=self.tensor3,
                                           k_grid=k_grid,
                                           sm1=1.0, sm0=1.0,
                                           nsm=1,
                                           q_path=self.qpath,
                                           T=T0)
# ----------
# Funciones
# ----------
def NombredeFuncion(arg):
    pass

def main(args):
    if (len( sys.argv ) > 1):
        #Fichero_dyn = "dyn_start_population25_"
        Fichero_dyn = arg[1]
        print ("Dyn_filename:",Fichero_dyn)
        #nqirr = 10
        nqirr = arg[2]
        print ("nqirr:",nqirr)
        #T0 = 50
        T0 = arg[3]
        print ("T0:",T0)
        #Espectro =  Funcion_espectral(Fichero_dyn,nqirr,PATH = "GXMGRX")
        Espectro =  Funcion_espectral(Fichero_dyn,nqirr)
        Espectro.prepara_tensor()
        starttime = timeit.default_timer()
        print("The start time is :",starttime)
        Espectro.calcula_espectro_basico_SrTiO3(T0)
        Espectro.calcula_espectro_basico_SrTiO3_multiprocessing(T0,4)
        Espectro.calcula_full_correction_en_punto_G(T0)
        Espectro.calcula_full_correction_en_punto_R(T0)
        Espectro.calcula_oneshot_correction_en_punto_Gamma(T0)
        Espectro.calcula_oneshot_correction_en_punto_R(T0)
        Espectro.calcula_oneshot_correction_along_PATH(T0)
        Espectro.calcula_espectro_correction_SrTiO3(T0)
        Espectro.calcula_espectro_correction_multiprocessing_SrTiO3(T0,8)
        print("The time difference is :", timeit.default_timer() - starttime)
        with open('Time_output.txt', 'w') as f:
            f.write("The time difference is :", timeit.default_timer() - starttime)
        Espectro.dibuja_espectro_basico_SrTiO3(filename = "SrTiO3_static.dat", PATH = "GXMGRX")
    else:
        print ("Arguments are dyn_filename, nqirr, T. In that order, separated by simple spaces")

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
