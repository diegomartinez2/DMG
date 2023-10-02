#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Lambda_test.py
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
#import Eliashberg
import Eliashberg_new as Eliashberg
import numpy as np
import os.path
import matplotlib.pyplot as plt
#from ase.units import create_units
from scipy import integrate

def main(arg):
    #units = create_units('2014')
    #file = 'texto2.txt'
    file = 'texto1.txt'
    if not (os.path.isfile(file)):
        print('Error, file not found')
    else:
        pars = np.loadtxt(file, usecols = (1,5,9,12))  #mode::Omega::Gamma::lambda

    superconductor = Eliashberg.Eliashberg_test(pars)
    Nef = 32.993055569433089 #np.loadtxt('DOS', usecols = (2)) #* factor_Nef
    #Nef = 6.060738
    #Nota en los datos de 'test' DOS = 6.060738 states/spin/Ry/Unit Cell
    print('Nef=',Nef,' ')
    print('lambda_q reded from data:',superconductor.pars[np.where(superconductor.pars==1),3])
    superconductor.energy = [0.0]
    #superconductor.Ne = [Nef]
    superconductor.Nef = Nef #* 0,0367493 from eV to Hartree #* superconductor.from_Ry_to_Hartree
    superconductor.pars[:,1] *= superconductor.from_cm1_to_Hartree
    superconductor.pars[:,2] *= superconductor.from_GHz_to_Hartree
    print("Omega range:",np.min(superconductor.pars[:,1]-np.abs(np.max(superconductor.pars[:,2]))),'::',np.max(superconductor.pars[:,1])+np.abs(np.max(superconductor.pars[:,2])))
    #frequencies = np.arange(0,20,1e-3) #    frequencies = np.arange(0,20,1e-7) <- it should have a plateau somewhere!! but with a step of 1e-7 still grows.
    frequencies = np.linspace(0,0.01,10000)
    lambda_1 = superconductor.Lambda_new(frequencies)
    print('Lambda_1=',lambda_1,':: test: must be 2.01') # Lambda calculated from Lambda_q
    #print('Lambda_2=',superconductor.lambda_2) #Lambda calculated fron Eliashberg function
    #print('Lambda_2_test',superconductor.lambda_2_test)
    print('Lambda_2_test2=',superconductor.lambda_2_test2,':: test: must be 2.03')
    np.savetxt('lambda_1_2.txt',(lambda_1,superconductor.lambda_2_test2))
    #print('Lambda_1/Lambda_2,Lambda_1/Lambda_2_test,Lambda_1/Lambda_2_test2:',lambda_1/superconductor.lambda_2,lambda_1/superconductor.lambda_2_test,lambda_1/superconductor.lambda_2_test2)
    #print('Lambda_1,Lambda_2*33,Lambda_2/33,Lambda_2_test*33,Lambda_2_test/33,Lambda_2_test2*33,Lambda_2_test2/33:',lambda_1,superconductor.lambda_2*33,superconductor.lambda_2/33,superconductor.lambda_2_test*33,superconductor.lambda_2_test/33,superconductor.lambda_2_test2*33,superconductor.lambda_2_test2/33)
#----------test--v----lambda_1---
    center = superconductor.pars[:,1] #* superconductor.from_cm1_to_Hartree
    width = superconductor.pars[:,2] #* superconductor.from_GHz_to_Hartree
    lamba = np.array([])
    for i in range(len(center)):
         lamba = np.append(lamba,superconductor.Lambda_q(width[i],center[i],Nef))
    np.savetxt('lambda_lista.txt',lamba)
    print('Test Gaussian->')
    suma = 0
    # for w in range(200000):
    #     suma += superconductor.test_gaussian(w/10000,superconductor.pars[:,1],len(superconductor.pars[:,1]))
    w = np.linspace(-100,100,20000)
    suma = integrate.simpson(superconductor.test_gaussian(w,superconductor.pars[:,1],len(superconductor.pars[:,1])), w)
    print(suma,':',suma/len(superconductor.pars[:,1]))
    pass



if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
