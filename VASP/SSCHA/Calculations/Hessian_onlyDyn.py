#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
#
#  Copyright 2024 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
from __future__ import print_function
from __future__ import division

# Import the modules to read the dynamical matrix
import cellconstructor as CC
import cellconstructor.Phonons

# Import the SCHA modules
import sscha, sscha.Ensemble


# Here the input information
DATA_DIR = 'pop3/data'               # path to the directory ens_pop#lastpop where the last population is stored
N_RANDOM = 512                                  # number elements in the ensamble
#DYN_PREFIX =  'pop3/dyn_start_population3_'          # dyn mat that generated the last population
#FINAL_DYN =   'pop3/dyn_end_population3_'            # SSCHA dyn mat obtained with the last minimization
DYN_file =   'pop3/dyn_end_population3_'            # SSCHA dyn mat obtained with the last minimization
SAVE_PREFIX = 'dyn_hessian_'                    # Free energy Hessian dynamical matrices output
NQIRR = 4                                       # The number or irredubcible q points
#Tg = 300                                         # The temperature used to generate the configurations
T =  300                                         # The temperature for the calculation
POPULATION = 3                                  # number of last population
INCLUDE_V4 = False                              # True to include the 4th-order SSCHA FC term to calculate the Hessian

# INFO = """
# In this example we compute the free energy hessian.
#
# The ensemble has been generated with the dynamical matrix at:
# {}
#
# And to compute the hessian we will use reweighting at:
# {}
#
# The original temperature was {} K, we use reweighting to {} K.
# The ensemble, population {}, is located at: {}
# The number of configuration is {}.
# Do we include the v4 in the calculation? {}
#
# The free energy Hessian will be saved in: {}
#
# The (symmetrized) 3rd order FCs in d3_realspace_sym.npy
#
# """.format(DYN_PREFIX, FINAL_DYN, Tg, T, POPULATION, DATA_DIR,
#            N_RANDOM, INCLUDE_V4, SAVE_PREFIX)


#print(INFO)
print()
print(" ======== RUNNING ======== ")
print()

# print("Loading the original dynamical matrix...")
# dyn = CC.Phonons.Phonons(DYN_PREFIX, NQIRR)
# print("Loading the current dynamical matrix...")
# final_dyn = CC.Phonons.Phonons(FINAL_DYN, NQIRR)
dyn = CC.Phonons.Phonons(DYN_file, NQIRR)

print("Loading the ensemble...")
#ens = sscha.Ensemble.Ensemble(dyn, Tg, dyn.GetSupercell())
ens = sscha.Ensemble.Ensemble(dyn, T, dyn.GetSupercell()) #See if we can use final_dyn here
ens.load(DATA_DIR, POPULATION, N_RANDOM)
# If the ensemble was saved in binary format, load it with
# ens.load_bin(DATA_DIR, POPULATION)

# print("Updating the importance sampling...")
# ens.update_weights(final_dyn, T)

print("Computing the free energy hessian...")
# Set get_full_hessian to false to have only the odd correction
# Usefull if you want to study the convergence with the number of configuration
dyn_hessian, d3  = ens.get_free_energy_hessian(include_v4 = INCLUDE_V4,
                                          get_full_hessian = True,
                                          verbose = True,
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
