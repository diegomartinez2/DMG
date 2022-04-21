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

ax1.set_xticks([0.0,0.707107,1.457107,1.707107,2.060660,2.560660,3.173033])
ax1.set_xticklabels(['X','$\Gamma$','K','M','W','L','$\Gamma$'])

ax1.set_xlim((0,3.173033))

ax1.grid(axis='x',color='black', linewidth=0.5)

ax1.set_ylim((0,125))

data = N.loadtxt('dyn_end_population3.freq.gp')
ax1.plot(data[:,0],data[:,1], color='r', label='pop3 auxiliary')
ax1.plot(data[:,0],data[:,2:7], color='r')

data = N.loadtxt('hessian.freq.gp')
ax1.plot(data[:,0],data[:,1], color='darkgreen', label='hessian')
ax1.plot(data[:,0],data[:,2:7], color='darkgreen')

ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)

ax1.axhspan(-40,0,color='lightgrey')

ax1.legend(frameon=False)

plt.savefig('hessian.pdf', bbox_inches='tight')

