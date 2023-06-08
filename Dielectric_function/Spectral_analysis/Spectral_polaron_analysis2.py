#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Spectral_polaron_analysis2.py
#
#  Copyright 2023 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
# ---------------------------
# Importación de los módulos
# ---------------------------
import matplotlib.pyplot as plt
import numpy as np
import scipy as scipy
from scipy.signal import find_peaks


class Polaron_analysis(object):
    """Code for the analysis of Silkin's polaron data."""

    def __init__(self, arg, namefile):
        super(Polaron_analysis, self).__init__()
        self.arg = arg
        self.namefile = namefile
        self.pars = np.zeros((51,3))
        self.perr_lorentz = np.zeros((51,3))

    def load_data(self):
        with open(self.namefile) as file:
             lines = [line.rsplit() for line in file]
        Spec = np.zeros(255051)
        for i in range(255051):
             Spec[i]=lines[i][6]
        #data = np.resize(Spec,(51,int(len(Spec)/51)))
        data = np.resize(Spec,(51,5001))
        Frequency = np.zeros(5001)
        for i in range(5001):
             Frequency[i]=lines[i][2]
        return data, Frequency

    def load_big_file(self, index, filename="A7_EPS.dat"):
        q_x, q_y, w, epsilon = np.loadtxt(filename,usecols=(0,1,2,6), unpack=True)
        print (len(np. unique(q_x))) #np. unique(my_array, return_counts=True)
        print (len(np. unique(q_y)))
        print (len(np. unique(w)))
        self.data = np.resize(epsilon,(51,51,5001))
        #plot.contour(data[i])
        """
        direccion q_y -> data[i]
        direccion q_x=q_y=i -> data = data[i][i] for i in {0..51}
        direccion
        """
        print (np.shape(self.data[int(index)]),"=(51,5001)?")
        return self.data[int(index)], w[:5001]

    def plot_contour(self,data):
        plt.style.use('_mpl-gallery-nogrid')
        fig, ax1 = plt.subplots(1,1)
        cax = ax1.imshow(data.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.3,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)

        cbar = fig.colorbar(cax)
        plt.show()
        pass

    def plot_contour2(self,data,data2):
        plt.style.use('_mpl-gallery-nogrid')
        fig, ax1 = plt.subplots(2,1)
        cax = ax1[0].imshow(data.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.3,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        cax2 = ax1[1].imshow(data2.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.3,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        #ax1[0].set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)
        #ax1[1].set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)
        print (self.pars[:,1])
        cax3 = ax1[1].scatter(x = range(len(self.pars[:,1])), y = self.pars[:,1]*5000,c='k',marker='x',s=10)
        cax3 = ax1[1].scatter(x = range(len(self.pars[:,1])), y = (self.pars[:,1]+self.pars[:,2])*5000,c='k',marker='1',s=10)
        cax3 = ax1[1].scatter(x = range(len(self.pars[:,1])), y = (self.pars[:,1]-self.pars[:,2])*5000,c='k',marker='2',s=10)

        cbar = fig.colorbar(cax)
        cbar2 = fig.colorbar(cax2)
        if (len( sys.argv) == 3):
            plt.savefig("Ajuste_{}".format("original_vs_fitted"))
        else:
            plt.savefig("Ajuste_{}_{}".format(self.namefile,"original_vs_fitted"))
        plt.show()
        pass

    def locate_1Lorenztian(self, Frequency, Spec, index):
        def _1Lorentzian(x, amp1, cen1, wid1):
            return amp1*wid1**2/((x-cen1)**2+wid1**2)
        peaks, _ = find_peaks(Spec,width=5,rel_height=0.3)
        print (peaks,Frequency[peaks])
        print (peaks,_["widths"])
        amp = 0.009
        cen = 0.009
        wid = 0.001
        popt_lorentz, pcov_lorentz = scipy.optimize.curve_fit(_1Lorentzian, Frequency, Spec, p0=[amp, cen, wid])

        perr_lorentz = np.sqrt(np.diag(pcov_lorentz))
        pars = popt_lorentz[:]
        lorentz_peak = _1Lorentzian(Frequency, *pars)
        print ("-------------Peak ",index,"-------------")
        print ("amplitude = %0.2f (+/-) %0.2f" % (pars[0], perr_lorentz[0]))
        print ("center = %0.2f (+/-) %0.2f" % (pars[1], perr_lorentz[1]))
        print ("width = %0.2f (+/-) %0.2f" % (pars[2], perr_lorentz[2]))
        print ("area = %0.2f" % np.trapz(lorentz_peak))
        print ("--------------------------------")

        # plt.plot(Frequency, Spec, label="original")
        # plt.plot(Frequency, lorentz_peak, "r",ls=':', label="ajuste")
        # #plt.plot(Frequency[peaks], Spec[peaks], "xr", label="peaks")
        # #plt.plot(Frequency[peaks], Spec[peaks], "1k", label="center")
        # plt.legend()
        # plt.xlim([0, Frequency.max()])
        # plt.ylim([0, Spec.max()])
        # plt.tight_layout()
        # plt.savefig("Ajuste_{}_{}".format(self.namefile,index))
        # plt.show()
        self.pars[index] = pars
        self.perr_lorentz[index] = perr_lorentz
        #np.savetxt(f, (pars[0], perr_lorentz[0],pars[1], perr_lorentz[1],pars[2], perr_lorentz[2]) , fmt='%1.3f', newline=", ")
        #f.write("\n")
        return lorentz_peak

    def fitting_Lorentz(self,frequencies,data):
        self.Fitted_data = np.zeros((51,5001))
        for i in range(51):
            self.Fitted_data[i] = self.locate_1Lorenztian(frequencies,data[i],i)
        pass

def main(arg):
    #namefile = "xbv"
    #print (len(sys.argv))
    if (len( sys.argv ) > 1):
        print ("nombre de fichero:",arg[1])
        namefile = arg[1]
        print ("Creando el objeto polaron")
        polaron = Polaron_analysis(arg,namefile)
        print ("Leyendo datos")

        if (len( sys.argv) == 3):
            index = arg[2]
            data, frequencies = polaron.load_big_file(index, arg[1])
        else:
            data, frequencies = polaron.load_data()

        print (np.shape(data),"=(51,5001)?")
        print ("dibuja")
        polaron.plot_contour(data)
        print("calcula ajuste loretzian")
        #none = polaron.locate_1Lorenztian(frequencies,data[30])
        #f=open('file_data_fittings.txt','a')
        #f.write("-amplitude--(+/-)error--center--(+/-)error--width--(+/-)error-\n")
        polaron.fitting_Lorentz(frequencies,data)
        #f.close()
        np.savetxt("data_fitting_amplitudes.txt",np.c_[polaron.pars[:,0], polaron.perr_lorentz[:,0]],header='#-----amplitude--(+/-)error---', footer='-------------')
        np.savetxt("data_fitting_center.txt",np.c_[polaron.pars[:,1], polaron.perr_lorentz[:,1]],header='#-----center(e.V.)--(+/-)error---', footer='-------------')
        np.savetxt("data_fitting_width.txt",np.c_[polaron.pars[:,2], polaron.perr_lorentz[:,2]],header='#-----width(e.V.)--(+/-)error---', footer='-------------')
        print ("dibuja")
        polaron.plot_contour(polaron.Fitted_data)
        polaron.plot_contour2(data,polaron.Fitted_data)
    else:
        print ("Arguments are namefile and the index of q_x as second argument if you want the BIG FILE")
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
