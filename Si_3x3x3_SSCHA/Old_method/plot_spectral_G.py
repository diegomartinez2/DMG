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

ax1.set_ylim((0,0.75))

data = N.loadtxt('nomm_spectral_func_1.00.dat')
ax1.plot(data[:,1],data[:,2], color='r', label='Full-no mode mixing')

data = N.loadtxt('nomm_spectral_func_lorentz_one_shot_1.00.dat')
ax1.plot(data[:,1],data[:,2], color='b', label='Lorentzian')


ax1.set_xlabel(r'Frequency (cm$^{-1}$)', fontsize=12)
ax1.set_xlabel(r'$\sigma(\omega)$)', fontsize=12)


ax1.legend(frameon=False)

plt.savefig('spectral_G.pdf', bbox_inches='tight')

