import Eliashberg
import numpy as np
import os.path
import matplotlib.pyplot as plt

def main(arg):
    #factor_Omega = 4.55633e-6 # from cm-1 to Hartree
    #factor_Gamma = 4.55633e-6/29.9793 # from GHz to Hartree
    #factor_Nef = 0.5 # from Ry to Hartree
    #file = 'texto2.txt'
    file = 'texto1.txt'
    if not (os.path.isfile(file)):
        print('Error, file not found')
    else:
        pars = np.loadtxt(file, usecols = (1,5,9))
    #print(pars)
    print('0--------------0')
    #print(pars[:,1])
    #print('Min(pars[0])=',np.amin(pars[:,1]))
    #pars[:,1] *= factor_Gamma
    #pars[:,2] *= factor_Omega

    superconductor = Eliashberg.Eliashberg(pars)
    Nef = np.loadtxt('DOS', usecols = (2)) #* factor_Nef
    print('Nef=',Nef)
    superconductor.energy = [0.0]
    superconductor.Ne = [Nef]
    superconductor.Nef = Nef * superconductor.from_Ry_to_Hartree
    superconductor.pars[:,1] *= superconductor.from_cm1_to_Hartree
    superconductor.pars[:,2] *= superconductor.from_GHz_to_Hartree
    #frequencies = np.arange(0,0.001,1e-9)
    frequencies = np.arange(0,1,1e-2) #    frequencies = np.arange(0,20,1e-7) <- it should have a plateau somewhere!! but with a step of 1e-7 still grows.
    lambda_1 = superconductor.Lambda(frequencies)
    print('Lambda_1=',lambda_1) # Lambda calculated from Lambda_q
    print('Lambda_2=',superconductor.lambda_2) #Lambda calculated fron Eliashberg function
    np.savetxt('lambda_1_2.txt',(lambda_1,superconductor.lambda_2))
#----------test--v----lambda_1---
    center = superconductor.pars[:,1] #* superconductor.from_cm1_to_Hartree
    width = superconductor.pars[:,2] #* superconductor.from_GHz_to_Hartree
    lamba = np.array([])
    for i in range(len(center)):
         lamba = np.append(lamba,superconductor.Lambda_q(width[i],center[i],Nef))
    np.savetxt('lambda_lista.txt',lamba)
    pass



if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))