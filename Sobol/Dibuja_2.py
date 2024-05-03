#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SSCHA_exercise.py
#
#  Copyright 2022 Diego <diego@u038025>
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

# Import the cellconstructor stuff
import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

# Import the modules of the force field
#import fforces as ff
#import fforces.Calculator

# Import the modules to run the sscha
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.Utilities

import spglib
from ase.visualize import view

# Import Matplotlib to plot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import timeit
def dibuja1():
    #----------------------SOBOL------------------------------------------------
        hessian_data = np.loadtxt("Sobol/16384_hessian_vs_temperature.dat")
        for i in range(len(hessian_data[:,2])):
            if (hessian_data[i,2]<0):
                hessian_data[i,2] = -(hessian_data[i,2]**2)
            else:
                hessian_data[i,2] = hessian_data[i,2]**2
        plt.plot(hessian_data[:,0], hessian_data[:,2], marker = "o")
        
    #----------------------RANDOM-----------------------------------------------
        hessian_data0 = np.loadtxt("Random/0_hessian_vs_temperature.dat")
        hessian_data1 = np.loadtxt("Random/1_hessian_vs_temperature.dat")
        hessian_data2 = np.loadtxt("Random/2_hessian_vs_temperature.dat")
        hessian_data3 = np.loadtxt("Random/3_hessian_vs_temperature.dat")
        hessian_data4 = np.loadtxt("Random/4_hessian_vs_temperature.dat")
        hessian_data5 = np.loadtxt("Random/5_hessian_vs_temperature.dat")

        a = [hessian_data0[:,2],hessian_data1[:,2],hessian_data2[:,2],hessian_data3[:,2],hessian_data4[:,2],hessian_data5[:,2]]
        media_data = sum(a)/6
        for i in range(len(media_data)):
            if (media_data[i]<0):
                media_data[i] = -(media_data[i]**2)
            else:
                media_data[i] = media_data[i]**2
        plt.plot(hessian_data2[:,0], media_data)
        max_a = np.zeros(len(hessian_data0[:,0]))
        min_a = np.zeros(len(hessian_data0[:,0]))
        for i in range(len(hessian_data0[:,0])):
            max_a[i]=max((a[0][i],a[1][i],a[2][i],a[3][i],a[4][i],a[5][i]))
            min_a[i]=min((a[0][i],a[1][i],a[2][i],a[3][i],a[4][i],a[5][i]))
        for i in range(len(max_a)):
            if (max_a[i]<0):
                max_a[i] = -(max_a[i]**2)
            else:
                max_a[i] = max_a[i]**2
            if (min_a[i]<0):
                min_a[i] = -(min_a[i]**2)
            else:
                min_a[i] = min_a[i]**2
        plt.fill_between(hessian_data2[:,0], max_a, min_a,alpha=0.2)
    #-------------------------SOBOL+SCATTER-------------------------------------
        hessian_data0 = np.loadtxt("SobolScatter/0_hessian_vs_temperature.dat")
        hessian_data1 = np.loadtxt("SobolScatter/1_hessian_vs_temperature.dat")
        hessian_data2 = np.loadtxt("SobolScatter/2_hessian_vs_temperature.dat")
        hessian_data3 = np.loadtxt("SobolScatter/3_hessian_vs_temperature.dat")
        hessian_data4 = np.loadtxt("SobolScatter/4_hessian_vs_temperature.dat")
        hessian_data5 = np.loadtxt("SobolScatter/5_hessian_vs_temperature.dat")

        a = [hessian_data0[:,2],hessian_data1[:,2],hessian_data2[:,2],hessian_data3[:,2],hessian_data4[:,2],hessian_data5[:,2]]
        media_data = sum(a)/6
        for i in range(len(media_data)):
            if (media_data[i]<0):
                media_data[i] = -(media_data[i]**2)
            else:
                media_data[i] = media_data[i]**2
        plt.plot(hessian_data2[:,0], media_data)
        max_a = np.zeros(len(hessian_data0[:,0]))
        min_a = np.zeros(len(hessian_data0[:,0]))
        for i in range(len(hessian_data0[:,0])):
            max_a[i]=max((a[0][i],a[1][i],a[2][i],a[3][i],a[4][i],a[5][i]))
            min_a[i]=min((a[0][i],a[1][i],a[2][i],a[3][i],a[4][i],a[5][i]))
        for i in range(len(max_a)):
            if (max_a[i]<0):
                max_a[i] = -(max_a[i]**2)
            else:
                max_a[i] = max_a[i]**2
            if (min_a[i]<0):
                min_a[i] = -(min_a[i]**2)
            else:
                min_a[i] = min_a[i]**2
        plt.fill_between(hessian_data2[:,0], max_a, min_a,alpha=0.2)

        #plt.figure(dpi = 120)
        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.xlabel("Temperature [K]")
        plt.ylabel("Frequency [cm-1]")
        plt.legend()
        plt.tight_layout()



        #plt.savefig('Temp_Freq_2.png')
        plt.show()

        # plt.figure(dpi = 120)
        # plt.plot(hessian_data[:,0], np.sign(hessian_data[:,2]) * hessian_data[:,2]**2, label = "Free energy curvature", marker = "o")
        # plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        # plt.xlabel("Temperature [K]")
        # plt.ylabel("$\omega^2$ [cm-2]")
        # plt.legend()
        # plt.tight_layout()
        # plt.savefig('Temp_Omeg.png'
        #plt.show()

def dibuja2():
        hessian_data = np.loadtxt("Sobol/16384_hessian_vs_configuraciones.dat")
        plt.plot(hessian_data[:,0], hessian_data[:,2], label = "Sobol Free energy curvature", marker = "o")
        hessian_data0 = np.loadtxt("Random/hessian_vs_configuraciones.dat")
        hessian_data1 = np.loadtxt("Random/1_hessian_vs_configuraciones.dat")
        hessian_data2 = np.loadtxt("Random/2_hessian_vs_configuraciones.dat")
        hessian_data3 = np.loadtxt("Random/3_hessian_vs_configuraciones.dat")
        hessian_data4 = np.loadtxt("Random/4_hessian_vs_configuraciones.dat")
        hessian_data5 = np.loadtxt("Random/5_hessian_vs_configuraciones.dat")
        #plt.figure(dpi = 120)
#        plt.plot(hessian_data[:,0], hessian_data[:,1], label = "Min SCHA freq", marker = ">")
#        plt.plot(hessian_data[:,0], hessian_data[:,2], label = "Free energy curvature", marker = "o")
#        plt.plot(hessian_data1[:,0], hessian_data1[:,2], marker = "o")
#        plt.plot(hessian_data2[:,0], hessian_data2[:,2], marker = "o")
#        plt.plot(hessian_data3[:,0], hessian_data3[:,2], marker = "o")
#        plt.plot(hessian_data4[:,0], hessian_data4[:,2], marker = "o")
#        plt.plot(hessian_data5[:,0], hessian_data5[:,2], marker = "o")
        a = [hessian_data0[:,2],hessian_data1[:,2],hessian_data2[:,2],hessian_data3[:,2],hessian_data4[:,2],hessian_data5[:,2]]
        media_data = sum(a)/6
        plt.plot(hessian_data2[:,0], media_data)
        max_a = np.zeros(len(hessian_data0[:,0]))
        min_a = np.zeros(len(hessian_data0[:,0]))
        for i in range(len(hessian_data0[:,0])):
            max_a[i]=max((a[0][i],a[1][i],a[2][i],a[3][i],a[4][i],a[5][i]))
            min_a[i]=min((a[0][i],a[1][i],a[2][i],a[3][i],a[4][i],a[5][i]))
        plt.fill_between(hessian_data2[:,0], max_a, min_a,alpha=0.2)

        hessian_data0 = np.loadtxt("0.001/hessian_vs_configuraciones.dat")
        hessian_data1 = np.loadtxt("0.001/1_hessian_vs_configuraciones.dat")
        hessian_data2 = np.loadtxt("0.001/2_hessian_vs_configuraciones.dat")
        hessian_data3 = np.loadtxt("0.001/3_hessian_vs_configuraciones.dat")
        hessian_data4 = np.loadtxt("0.001/4_hessian_vs_configuraciones.dat")
        hessian_data5 = np.loadtxt("0.001/5_hessian_vs_configuraciones.dat")
        #plt.figure(dpi = 120)
#        plt.plot(hessian_data[:,0], hessian_data[:,1], label = "Min SCHA freq", marker = ">")
#        plt.plot(hessian_data[:,0], hessian_data[:,2], label = "Free energy curvature", marker = "o")
#        plt.plot(hessian_data1[:,0], hessian_data1[:,2], marker = "o")
#        plt.plot(hessian_data2[:,0], hessian_data2[:,2], marker = "o")
#        plt.plot(hessian_data3[:,0], hessian_data3[:,2], marker = "o")
#        plt.plot(hessian_data4[:,0], hessian_data4[:,2], marker = "o")
#        plt.plot(hessian_data5[:,0], hessian_data5[:,2], marker = "o")
        a = [hessian_data0[:,2],hessian_data1[:,2],hessian_data2[:,2],hessian_data3[:,2],hessian_data4[:,2],hessian_data5[:,2]]
        media_data = sum(a)/6
        plt.plot(hessian_data2[:,0], media_data)
        max_a = np.zeros(len(hessian_data0[:,0]))
        min_a = np.zeros(len(hessian_data0[:,0]))
        for i in range(len(hessian_data0[:,0])):
            max_a[i]=max((a[0][i],a[1][i],a[2][i],a[3][i],a[4][i],a[5][i]))
            min_a[i]=min((a[0][i],a[1][i],a[2][i],a[3][i],a[4][i],a[5][i]))
        plt.fill_between(hessian_data2[:,0], max_a, min_a,alpha=0.2)

        plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        plt.axvline(x=1024, ls = "dotted")
        plt.axvline(x=512, ls = "dotted")
        plt.axvline(x=256, ls = "dotted")
        plt.axvline(x=128, ls = "dotted")
        plt.axvline(x=64, ls = "dotted")
        plt.axvline(x=32, ls = "dotted")
        plt.xlabel("Num. Configurations")
        plt.ylabel("Frequency [cm-1]")
        plt.legend()
        plt.tight_layout()
        plt.savefig('Conf_Freq_2.png')
        plt.show()

        # plt.figure(dpi = 120)
        # plt.plot(hessian_data[:,0], np.sign(hessian_data[:,2]) * hessian_data[:,2]**2, label = "Free energy curvature", marker = "o")
        # plt.axhline(0, 0, 1, color = "k", ls = "dotted") # Draw the zero
        # plt.xlabel("Num. Configurations")
        # plt.ylabel("$\omega^2$ [cm-2]")
        # plt.legend()
        # plt.tight_layout()
        # plt.savefig('Conf_Omeg.png')
        #plt.show()

def estimate_coef(x, y):
	# number of observations/points
	n = np.size(x)

	# mean of x and y vector
	m_x = np.mean(x)
	m_y = np.mean(y)

	# calculating cross-deviation and deviation about x
	SS_xy = np.sum(y*x) - n*m_y*m_x
	SS_xx = np.sum(x*x) - n*m_x*m_x

	# calculating regression coefficients
	b_1 = SS_xy / SS_xx
	b_0 = m_y - b_1*m_x

	return (b_0, b_1)

def plot_regression_line(x, y, b):
	# plotting the actual points as scatter plot
	plt.scatter(x, y, color = "m",
			marker = "o", s = 30)

	# predicted response vector
	y_pred = b[0] + b[1]*x

	# plotting the regression line
	plt.plot(x, y_pred, color = "g")

	# plot line at y=0
	plt.axhline(y=0.0, color='r', linestyle='--')

	# putting labels
	plt.xlabel('T(K)')
	plt.ylabel('$\Omega^2$ (cm-1)')
	plt.title(r'Squared Frequencies versus Temperatures', fontsize='small')
	plt.tight_layout()
	plt.savefig("Linear_regression.png")
	#plt.savefig("Ajuste_{}".format("lineal"))
	# function to show plot
	plt.show()


def main(args):
    dibuja1()
    dibuja2()
    # b = estimate_coef(x, y)
	# print("Estimated coefficients:\nb_0 = {} \
	# 	\nb_1 = {}".format(b[0], b[1]))
    #
	# # plotting regression line
	# plot_regression_line(x, y, b)
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
