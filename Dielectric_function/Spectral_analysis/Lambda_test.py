import Spectral_plasmons_analysis
import numpy as np
import os.path
import matplotlib.pyplot as plt

def main(arg):
    factor_Omega = 4.55633e-6 # from cm-1 to Hartree
    factor_Gamma = 4.55633e-6/29.9793 # from GHz to Hartree
    file = 'texto2.txt'
    if not (os.path.isfile(file)):
        print('Error, file not found')
    else:
        pars = np.loadtxt(file, usecols = (1,5,9))
    print(pars)
    print('0--------------0')
    print(pars[:,1])
    print('Min(pars[0])=',np.amin(pars[:,1]))
    pars[:,1] *= factor_Gamma
    pars[:,2] *= factor_Omega

    superconductor = Spectral_plasmons_analysis.Eliashberg(pars)
    Nef = np.loadtxt('DOS', usecols = (2))
    print('Nef=',Nef)
    superconductor.energy = [0.0]
    superconductor.Ne = [Nef]
    superconductor.Nef = Nef
    frequencies = np.arange(0,0.01,0.0001)
    lambda_1 = superconductor.Lambda(frequencies)
    print('Lambda_1=',lambda_1) # Lambda calculated from Lambda_q
    print('Lambda_2=',superconductor.lambda_2) #Lambda calculated fron Eliashberg function
    pass



if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
