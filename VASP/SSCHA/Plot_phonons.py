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
ax1.set_ylim((-40,125))

ax1.grid(axis='x',color='black', linewidth=0.5)

data = N.loadtxt('harmonic_calculation/harmonic.freq.gp')
ax1.plot(data[:,0],data[:,1], color='k', linestyle='dotted', label='harmonic')
ax1.plot(data[:,0],data[:,2:7], color='k', linestyle='dotted')

data = N.loadtxt('dyn_start_population1.freq.gp')
ax1.plot(data[:,0],data[:,1], color='k', linestyle='dashed', label='starting')
ax1.plot(data[:,0],data[:,2:7], color='k', linestyle='dashed')

data = N.loadtxt('dyn_end_population1.freq.gp')
ax1.plot(data[:,0],data[:,1], color='b', label='pop1')
ax1.plot(data[:,0],data[:,2:7], color='b')

data = N.loadtxt('dyn_end_population2.freq.gp')
ax1.plot(data[:,0],data[:,1], color='g', label='pop2')
ax1.plot(data[:,0],data[:,2:7], color='g')

data = N.loadtxt('dyn_end_population3.freq.gp')
ax1.plot(data[:,0],data[:,1], color='r', label='pop3')
ax1.plot(data[:,0],data[:,2:7], color='r')
for i in range(1,3,1):
    data = N.loadtxt('dyn_end_population{i}.freq.gp')
    ax1.plot(data[:,0],data[:,1], color='r', label='pop{i}'
    ax1.plot(data[:,0],data[:,2:7], color='r')

ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)

ax1.axhspan(-40,0,color='lightgrey')

ax1.legend(frameon=False)

plt.savefig('sscha_phonon_evolution.pdf', bbox_inches='tight')
