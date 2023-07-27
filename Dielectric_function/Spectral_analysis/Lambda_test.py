import Spectral_plasmons_analysis
import numpy as np
import os.path
import matplotlib.pyplot as plt

def main(arg):
    file = 'texto2.txt'
    if not (os.path.isfile(file)):
        print('Error, file not found')
    else:
        pars = np.loadtxt(file, usecols = (1,5,9))
    print(pars)
    print('0--------------0')
    print(pars[:,0])
    print('Min(pars[0])=',np.amin(pars[:,0]))
    superconductor = Spectral_plasmons_analysis.Eliashberg(pars)
    Nef = np.loadtxt('DOS', usecols = (2))
    print('Nef=',Nef)
    superconductor.energy = [0.0]
    superconductor.Ne = [Nef]
    superconductor.Nef = Nef
    frequencies = np.arange(0,1000,1)
    lambda_1 = superconductor.Lambda(frequencies)
    print('Lambda_1=',lambda_1) # Lambda calculated from Lambda_q
    print('Lambda_2=',superconductor.lambda_2) #Lambda calculated fron Eliashberg function
    pass



if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
