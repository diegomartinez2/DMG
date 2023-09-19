import Eliashberg
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
        pars = np.loadtxt(file, usecols = (1,5,9))

    superconductor = Eliashberg.Eliashberg(pars)
    Nef = 32.993055569433089 #np.loadtxt('DOS', usecols = (2)) #* factor_Nef
    print('Nef=',Nef,' ?=33?')
    superconductor.energy = [0.0]
    superconductor.Ne = [Nef]
    superconductor.Nef = Nef *0.367493# 0,0367493 from eV to Hartree #* superconductor.from_Ry_to_Hartree
    superconductor.pars[:,1] *= superconductor.from_cm1_to_Hartree
    superconductor.pars[:,2] *= superconductor.from_GHz_to_Hartree
    # superconductor.pars[:,1] *= 1 #if using Hartree,
    # superconductor.pars[:,2] *= 1
    #frequencies = np.arange(0,20,1e-3) #    frequencies = np.arange(0,20,1e-7) <- it should have a plateau somewhere!! but with a step of 1e-7 still grows.
    frequencies = np.linspace(-20,20,200000)
    lambda_1 = superconductor.Lambda_new(frequencies)
    lambda_1 *= Nef #misterious factor... joking, this is the number of nodes in the example.
    print('Lambda_1=',lambda_1) # Lambda calculated from Lambda_q
    #print('Lambda_2=',superconductor.lambda_2) #Lambda calculated fron Eliashberg function
    print('Lambda_2_test',superconductor.lambda_2_test)
    print('Lambda_2_test2',superconductor.lambda_2_test2)
    np.savetxt('lambda_1_2.txt',(lambda_1,superconductor.lambda_2))
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
