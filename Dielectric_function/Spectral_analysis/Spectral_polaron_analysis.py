#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
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
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy as scipy
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import matplotlib.ticker as ticker
# -------
# Clases
# -------
class Polaron_analysis(object):
    """docstring for NombredeClase."""

    def __init__(self, arg):
        #super(NombredeClase, self).__init__()
        self.arg = arg

        self.Spec = np.zeros((51,51,5001))


    def read_data(self):
        with open("A7_EPS.dat") as file:
             lines = [line.rsplit() for line in file]
        for i in range(51*51*5001):
                    ii=int(lines[i][0]/0.00869)
                    ij=int(lines[i][1]/0.00869)
                    ik=int(lines[i][2]/0.0002)
                    print (i,ii,ij,ik)
                    self.Spec[ii][ij][ik]=lines[index][6]
        # for i in range(5001*50*51):
        #     if lines[i][0]=0.31290
        #     for j in range(50):
        #
        #      Spec[j][i]=lines[i][6]
        return 0

    def drawn_3D(self, X, Y, Z):
        levels = np.linspace(Z.min(), Z.max(), 7)

        # plot
        fig, ax = plt.subplots()

        ax.contourf(X, Y, Z, levels=levels)

        plt.show()
        return 0

    def drawn_3D2(self, data):
        #data = np.resize(Spec,(51,int(len(Spec)/51)))
        fig, ax1 = plt.subplots(1,1)
        cax = ax1.imshow(data.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.3,
	              cmap=plt.colormaps['jet'], origin='lower',
	              interpolation='gaussian', aspect='auto')
        ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)

        cbar = fig.colorbar(cax)
        plt.show()
        return 0

    def drawn_2D(self, Frequency, Spec, Data5):
        plt.plot(Frequency, Spec, label="original")
        plt.plot(Frequency, Data5*0.000006, label="Im(Epsilon)", ls=":")
        plt.plot(Frequency, Spec-(Data5*0.000006), label="Substraction")
        plt.xlim([0, Frequency.max()])
        plt.ylim([0, Spec.max()])
        plt.legend()
        plt.tight_layout()
        #plt.savefig("Epsilon_{}".format(namefile))
        plt.show()
        return 0

    def remove_noise(self):
        return 0

    def locate_1Lorenztian(self, Frequency, Spec):
        def _1Lorentzian(x, amp1, cen1, wid1):
            return amp1*wid1**2/((x-cen1)**2+wid1**2)
        popt_lorentz, pcov_lorentz = scipy.optimize.curve_fit(_1Lorentzian, Frequency, Spec, p0=[amp, cen, wid])

        perr_3lorentz = np.sqrt(np.diag(pcov_lorentz))
        peaks, _ = find_peaks(Spec,width=5,rel_height=0.3)
        print (peaks,Frequency[peaks])
        print (peaks,_["widths"])
        return 0

    def locate_3Lorenztians(self,amp1, cen1, wid1,amp2, cen2, wid2, amp3, cen3, wid3):
        def _1Lorentzian(x, amp, cen, wid):
            return amp*wid**2/((x-cen)**2+wid**2)
        def _3Lorentzian(x, amp1, cen1, wid1, amp2,cen2,wid2, amp3,cen3,wid3):
            return (amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
                    (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                        (amp3*wid3**2/((x-cen3)**2+wid3**2))

        popt_3lorentz, pcov_3lorentz = scipy.optimize.curve_fit(_3Lorentzian, x_array, y_array_3lorentz, p0=[amp1, cen1, wid1, \
                                                                                            amp2, cen2, wid2, amp3, cen3, wid3])

        perr_3lorentz = np.sqrt(np.diag(pcov_3lorentz))

        pars_1 = popt_3lorentz[0:3]
        pars_2 = popt_3lorentz[3:6]
        pars_3 = popt_3lorentz[6:9]
        lorentz_peak_1 = _1Lorentzian(x_array, *pars_1)
        lorentz_peak_2 = _1Lorentzian(x_array, *pars_2)
        lorentz_peak_3 = _1Lorentzian(x_array, *pars_3)

        # this cell prints the fitting parameters with their errors
        print ("-------------Peak 1-------------")
        print ("amplitude = %0.2f (+/-) %0.2f" % (pars_1[0], perr_3lorentz[0]))
        print ("center = %0.2f (+/-) %0.2f" % (pars_1[1], perr_3lorentz[1]))
        print ("width = %0.2f (+/-) %0.2f" % (pars_1[2], perr_3lorentz[2]))
        print ("area = %0.2f" % np.trapz(lorentz_peak_1))
        print ("--------------------------------")
        print ("-------------Peak 2-------------")
        print ("amplitude = %0.2f (+/-) %0.2f" % (pars_2[0], perr_3lorentz[3]))
        print ("center = %0.2f (+/-) %0.2f" % (pars_2[1], perr_3lorentz[4]))
        print ("width = %0.2f (+/-) %0.2f" % (pars_2[2], perr_3lorentz[5]))
        print ("area = %0.2f" % np.trapz(lorentz_peak_2))
        print ("--------------------------------")
        print ("-------------Peak 3-------------")
        print ("amplitude = %0.2f (+/-) %0.2f" % (pars_3[0], perr_3lorentz[6]))
        print ("center = %0.2f (+/-) %0.2f" % (pars_3[1], perr_3lorentz[7]))
        print ("width = %0.2f (+/-) %0.2f" % (pars_3[2], perr_3lorentz[8]))
        print ("area = %0.2f" % np.trapz(lorentz_peak_3))
        print ("--------------------------------")

        residual_3lorentz = y_array_3lorentz - (_3Lorentzian(x_array, *popt_3lorentz))
        return 0

    def  locate_1_to_3Lorenztians(self,amp1, cen1, wid1,amp2, cen2, wid2, amp3, cen3, wid3):
        peaks, _ = find_peaks(Spec, prominence=1)
        m = np.zeros(Frequency.shape, dtype=bool)
        m[peaks] = True
        x_max_range=Frequency.max()
        print (x_max_range)
        #--------------------------------
        x_array = Frequency
        y_array_3lorentz = Spec
        #--------------------------------
        def _base(x):
            if (x < 0.0):
        #    if (x < 0.0271):
                x_out = 0
            elif (x <= 0.0271):
                x_out = -1.9985*x*x+0.0651067*x+0.000119668
            elif ((x >= 0.0271)and(x <= 0.107)):
                x_out = -1.99329*x*x+0.264442*x-0.00527247
            else:
                x_out = 0
            return x_out

        def _1Lorentzian(x, amp1, cen1, wid1):
            return amp1*wid1**2/((x-cen1)**2+wid1**2)

        def _3Lorentzian(x, amp1, cen1, wid1, amp2,cen2,wid2, amp3,cen3,wid3):
                    return (amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
                            (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                                (amp3*wid3**2/((x-cen3)**2+wid3**2))

        def _4Lorentzian(x, amp1, cen1, wid1, amp2,cen2,wid2, amp3,cen3,wid3, amp4,cen4,wid4):
                    return (amp1*wid1**2/((x-cen1)**2+wid1**2)) + (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                                (amp3*wid3**2/((x-cen3)**2+wid3**2)) + (amp4*wid4**2/((x-cen4)**2+wid4**2))

        #--------------------------------------------
        plt.plot(Frequency, Spec, label="original")
        base = np.zeros(Frequency.shape)
        for i in range(len(Frequency)):
           base[i]=_base(Frequency[i])
        plt.plot(Frequency,base, label="base")
        plt.xlim([0, x_max_range])
        plt.ylim([0,0.015])
        plt.tight_layout()
        plt.savefig("Original_{}".format(namefile))
        plt.show()
        #--------------------------------------------
        if (amp2==0.0):
            popt_1lorentz, pcov_1lorentz = scipy.optimize.curve_fit(_1Lorentzian, x_array, (y_array_3lorentz - base), p0=[amp1, cen1, wid1], bounds=(0,np.inf))
            pars_1 = popt_1lorentz[:]
            perr_1lorentz = np.sqrt(np.diag(pcov_1lorentz))

        else:
            popt_3lorentz, pcov_3lorentz = scipy.optimize.curve_fit(_3Lorentzian, x_array, (y_array_3lorentz - base), p0=[amp1, cen1, wid1, amp2, cen2, wid2, amp3, cen3, wid3], bounds=(0,np.inf))
            pars_1 = popt_3lorentz[0:3]
            pars_2 = popt_3lorentz[3:6]
            pars_3 = popt_3lorentz[6:9]
            perr_3lorentz = np.sqrt(np.diag(pcov_3lorentz))

        #popt_4lorentz, pcov_4lorentz = scipy.optimize.curve_fit(_4Lorentzian, x_array, y_array_3lorentz, p0=[amp1, cen1, wid1, amp2, cen2, wid2, amp3, cen3, wid3, amp4, cen4, wid4])

        #perr_4lorentz = np.sqrt(np.diag(pcov_4lorentz))


        #pars_1 = popt_4lorentz[0:3]
        #pars_2 = popt_4lorentz[3:6]
        #pars_3 = popt_4lorentz[6:9]
        #pars_4 = popt_4lorentz[9:12]
        lorentz_peak_1 = _1Lorentzian(x_array, *pars_1)
        if (amp2!=0.0):
            lorentz_peak_2 = _1Lorentzian(x_array, *pars_2)
            lorentz_peak_3 = _1Lorentzian(x_array, *pars_3)
        #lorentz_peak_4 = _1Lorentzian(x_array, *pars_4)

        ## this cell prints the fitting parameters with their errors
        if (amp2==0.0):
            print ("-------------Peak 1-------------")
            print ("amplitude = {} (+/-) {}".format(pars_1[0], perr_1lorentz[0]))
            print ("center = {} (+/-) {}".format(pars_1[1], perr_1lorentz[1]))
            print ("width = {} (+/-) {}".format(pars_1[2], perr_1lorentz[2]))
            print ("area = {}".format(np.trapz(lorentz_peak_1)))
            print ("--------------------------------")
        else:
            print ("-------------Peak 1-------------")
            print ("amplitude = {} (+/-) {}".format(pars_1[0], perr_3lorentz[0]))
            print ("center = {} (+/-) {}".format(pars_1[1], perr_3lorentz[1]))
            print ("width = {} (+/-) {}".format(pars_1[2], perr_3lorentz[2]))
            print ("area = {}".format(np.trapz(lorentz_peak_1)))
            print ("--------------------------------")
            print ("-------------Peak 2-------------")
            print ("amplitude = {} (+/-) {}".format(pars_2[0], perr_3lorentz[3]))
            print ("center = {} (+/-) {}".format(pars_2[1], perr_3lorentz[4]))
            print ("width = {} (+/-) {}".format(pars_2[2], perr_3lorentz[5]))
            print ("area = {}".format(np.trapz(lorentz_peak_2)))
            print ("--------------------------------")
            print ("-------------Peak 3-------------")
            print ("amplitude = {} (+/-) {}".format(pars_3[0], perr_3lorentz[6]))
            print ("center = {} (+/-) {}".format(pars_3[1], perr_3lorentz[7]))
            print ("width = {} (+/-) {}".format(pars_3[2], perr_3lorentz[8]))
            print ("area = {}".format(np.trapz(lorentz_peak_3)))
            print ("--------------------------------")

        #print ("-------------Peak 1-------------")
        #print ("amplitude = %0.2f (+/-) %0.2f" % (pars_1[0], perr_4lorentz[0]))
        #print ("center = %0.2f (+/-) %0.2f" % (pars_1[1], perr_4lorentz[1]))
        #print ("width = %0.2f (+/-) %0.2f" % (pars_1[2], perr_4lorentz[2]))
        #print ("area = %0.2f" % np.trapz(lorentz_peak_1))
        #print ("--------------------------------")
        #print ("-------------Peak 2-------------")
        #print ("amplitude = %0.2f (+/-) %0.2f" % (pars_2[0], perr_4lorentz[3]))
        #print ("center = %0.2f (+/-) %0.2f" % (pars_2[1], perr_4lorentz[4]))
        #print ("width = %0.2f (+/-) %0.2f" % (pars_2[2], perr_4lorentz[5]))
        #print ("area = %0.2f" % np.trapz(lorentz_peak_2))
        #print ("--------------------------------")
        #print ("-------------Peak 3-------------")
        #print ("amplitude = %0.2f (+/-) %0.2f" % (pars_3[0], perr_4lorentz[6]))
        #print ("center = %0.2f (+/-) %0.2f" % (pars_3[1], perr_4lorentz[7]))
        #print ("width = %0.2f (+/-) %0.2f" % (pars_3[2], perr_4lorentz[8]))
        #print ("area = %0.2f" % np.trapz(lorentz_peak_3))
        #print ("--------------------------------")
        #print ("-------------Peak 4-------------")
        #print ("amplitude = %0.2f (+/-) %0.2f" % (pars_4[0], perr_4lorentz[9]))
        #print ("center = %0.2f (+/-) %0.2f" % (pars_4[1], perr_4lorentz[10]))
        #print ("width = %0.2f (+/-) %0.2f" % (pars_4[2], perr_4lorentz[11]))
        #print ("area = %0.2f" % np.trapz(lorentz_peak_4))
        #print ("--------------------------------")
        if (amp2==0.0):
            residual_1lorentz = y_array_3lorentz - (_1Lorentzian(x_array, *popt_1lorentz))
        else:
            residual_3lorentz = y_array_3lorentz - (_3Lorentzian(x_array, *popt_3lorentz))
        #residual_4lorentz = y_array_3lorentz - (_4Lorentzian(x_array, *popt_4lorentz))
        #------------------------------


        plt.plot(Frequency, Spec, label="original")
        if (amp2==0.0):
            plt.plot(Frequency, lorentz_peak_1+base, "r",ls=':', label="ajuste")
        else:
            plt.plot(Frequency, lorentz_peak_1+lorentz_peak_2+lorentz_peak_3+base, "r",ls=':', label="ajuste")
        #plt.plot(Frequency, lorentz_peak_1+lorentz_peak_2+lorentz_peak_3+lorentz_peak_4, "ro", label="ajuste")
        plt.plot(Frequency, base, label="base")
        plt.plot(Frequency, lorentz_peak_1, label="1")
        if (amp2!=0):
            plt.plot(Frequency, lorentz_peak_2, label="2")
            plt.plot(Frequency, lorentz_peak_3, label="3")
        #plt.plot(Frequency, lorentz_peak_4, label="4")
        plt.legend()
        plt.xlim([0, x_max_range])
        plt.ylim([0,0.015])
        plt.tight_layout()
        plt.savefig("Ajuste_{}".format(namefile))
        plt.show()

class Plot_3D(object):
    def __init__(self):
        return 0


# ----------
# Funciones
# ----------
def Detecta_picos(arg):
    #from peakdetect import peakdetect

    peaks = peakdetect(data, lookahead=20)
    # Lookahead is the distance to look ahead from a peak to determine if it is the actual peak.
    # Change lookahead as necessary
    higherPeaks = np.array(peaks[0])
    lowerPeaks = np.array(peaks[1])
    plt.plot(data)
    plt.plot(higherPeaks[:,0], higherPeaks[:,1], 'ro')
    plt.plot(lowerPeaks[:,0], lowerPeaks[:,1], 'ko')
    pass

def main(args):
    Spec = np.zeros(5001)
    Data5 = np.zeros(5001)
    polaron = Polaron_analysis(args)
    polaron.read_data()
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
