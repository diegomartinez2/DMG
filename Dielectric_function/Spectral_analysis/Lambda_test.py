import Spectral_plasmons_analysis
import numpy as np
import os.path
import matplotlib.pyplot as plt

def main(arg):
    file = 'text.txt'
    if not (os.path.isfile(file)):
        print('Error, file not found')
    else:
        pars = np.loadtxt(file, usecols = (6,10,13))
    pass



if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
