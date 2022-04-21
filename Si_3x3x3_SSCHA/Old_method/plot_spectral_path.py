#!/usr/local/bin/python

import numpy as N
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib import cm
import subprocess as SP
from shutil import copyfile

# Prepare plot of phonon spectra

fig, ax1 = plt.subplots(1,1)

data = N.loadtxt('nomm_spectral_func_in_path_1.00.dat')

plt.scatter(data[:,0], data[:,1], c=data[:,2], cmap='hot')

ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)

plt.colorbar()

plt.savefig('spectral_path.pdf', bbox_inches='tight')

