#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Spectral_plasmons_analysis.py
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


class Plasmon_analysis(object):
    """Code for the analysis of Silkin's plasmon data."""

    def __init__(self, arg, namefile):
        super(Plasmon_analysis, self).__init__()
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
        print ("q_x=",lines[0][0])
        return data, Frequency, lines[0][0]

    def load_big_file(self, index, filename="A7_EPS.dat", diagonal=False):
        q_x, q_y, w, epsilon = np.loadtxt(filename,usecols=(0,1,2,6), unpack=True)
        print (len(np. unique(q_x))) #np. unique(my_array, return_counts=True)
        print (len(np. unique(q_y)))
        print (len(np. unique(w)))
        #self.data = np.resize(epsilon,(51,51,5001))
        self.data = np.resize(epsilon,(len(np. unique(q_x)),
                                        len(np. unique(q_y)),
                                        len(np. unique(w))))
        #plot.contour(data[i])
        """
        direccion q_y -> data[i]
        direccion q_x=q_y=i -> data = data[i][i] for i in {0..51}
        direccion
        """
        data = np.zeros((len(np. unique(q_x)),len(np. unique(w))))
        print ("+++++")
        if (diagonal == True):
            for i in range(len(np. unique(q_x))):
                print ("[" + "=" * i +  " " * ((len(np. unique(q_x)) - i)) + "]" +  str(i) + "%",end='\r', flush=True)
                for j in range(len(np. unique(w))):
                    data[i,j] = self.data[i,i,j]

        print ("*************")
        print (data.shape)
        # print (np.diagonal(self.data,axis1=0, axis2=1).T)
        # print (np.diagonal(self.data,axis1=0, axis2=1).shape)
        print ("--------------")

        print (np.shape(self.data[int(index)]),"=(51,5001)?")
        if (diagonal == True):
            return data, w[:5001]
        else:
            return self.data[int(index)], w[:5001]

    def plot_contour(self,data):
        plt.style.use('_mpl-gallery-nogrid')
        fig, ax1 = plt.subplots(1,1)
        fig.set_size_inches(10, 5)
        fig.set_dpi(100)
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
        fig.set_size_inches(10, 5)
        fig.set_dpi(100)
        cax = ax1[0].imshow(data.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.175,
        	vmin = 0.0 , vmax = 0.12,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        cax2 = ax1[1].imshow(data2.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.3,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
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

class Eliashberg(object):
    """docstring for Eliashberg."""

    def __init__(self, pars):
        #super(Eliashberg, self).__init__()
        self.pars = pars

    def e_spectral_function(self):
        print ("shape = ",np.shape(self.pars))
        center = self.pars[:,1]
        width = self.pars[:,2]
        N = 2 #electron spin up and down
        factor1 = 1/(2*np.pi*self.Ne*N)
        for i in range(len(center)):
            suma += (width/center) * self.gaussian(center)
        return factor1*summa

    def read_Ne(self,filename="out_DOS.dat"):
        self.energy, self.Ne = np.loadtxt(filename,usecols=(0,1), unpack=True)
        #self.Ne_0[np.where(energy==0.0)]
        numero_de_elementos = np.trapz(self.Ne[:(np.where(self.energy==0.0)[0][0])+1],dx=(self.energy[9]-self.energy[0])/10)
        print ("Density of states at Fermi level per cell=",self.Ne[np.where(self.energy==0.0)])
        print ("Number of elements (for the factor)=",numero_de_elementos)
        pass

    def Lambda_q(self,gamma_q,omega_q,Nef):
        Lamb_q=(1/(np.pi*Nef)) * (gamma_q/omega_q**2) #fix from omega to omega²
        return Lamb_q

    def a2F(self):
        print ("shape = ",np.shape(self.pars))
        center = self.pars[:,1]
        width = self.pars[:,2]
        #factor1 = 1/(2*self.Ne[np.where(self.energy==0.0)])
        factor1 = 1/len(center)     #test
        summa = 0
        for i in range(len(center)):
            summa += (width*center) * self.gaussian(x,center)
            # a2F(i) + gaussian * lambda(j,l) * w(j,l) * weightq(j) &
            #           / (2.0d0 * dble(total_qpoints_a2f))
        return factor1*summa # /2 or *0.5 ???

    def gaussian(self,x, center,gauss_width):
        #d = self.pars[:,0] #amplitude
        #center = self.pars[:,1]
        mu = center #self.pars[:1] #center
        sigma = gauss_width #self.pars[:,2] #width
        d = 1/(np.sqrt(2*np.pi)*sigma)
        g = d*np.exp(-( (x-mu)**2 / ( 2.0 * sigma**2 ) ) )
        #exp(-(w_aux-w(j,l))**2.0d0/(2.0d0*broad**2.0d0))/ (broad*sqrt(twopi)) #copy from qe-5.1.0_elph
        return g

    def Lambda(self,Frequencies):
        """
        Calculates the Lambda by two methods, notice that it must calculate the integral
        in a range that takes the Lorenztian obtained by the Plasmon_analysis object.
        """
        center = self.pars[:,1]
        width = self.pars[:,2]
        Ne = self.Ne[np.where(self.energy==0.0)]
        center = np.absolute(center) #test to force the abs
        width = np.absolute(width)
        summa1 = 0
        for i in range(len(center)):
            summa1 += self.Lambda_q(width,center,Ne)
        #Lambda_1=summa1/Ne
        Lambda_1=summa1/len(center)
        self.lambda_2=[]
        Frequncies = Frequencies+1
        Frequencies = np.append(Frequencies,Frequncies, axis=0)
        # for x in Frequencies:
        #     if (x == 0):
        #         continue
        #     else:
        #         a2F_x = self.a2F_2(x)/x
        #     print('Integrando(',x,')')
        #     self.lambda_2.append(2*np.trapz(a2F_x,dx=(Frequencies[9]-Frequencies[0])/10))  #<--This is OK?? (I think not, the integral is on frequencies -outside of the loop-)
        #     #self.lambda_2.append(2*np.trapz(a2F_x,Frequencies[1:])) #test this way of doing the integral
        #     #self.lambda_2.append(2*np.trapz(a2F_x)) #test this way of doing the integral <-NO
        for w in Frequencies:
            print('Frequence(',w,')              ', end="\r", flush=True)
            if (w == 0):
                a2F_x = []
            else:
                a2F_x.append(self.a2F_2(w)/w)
        self.lambda_2 = 2*np.trapz(a2F_x,dx=(Frequencies[9]-Frequencies[0])/10)
        print("len(freq[1:]),len(lambda1),len(lambda2),len(a2F_x) and shape of a2F_x")
        print(len(Frequencies[1:]),len(Lambda_1),len(self.lambda_2),len(a2F_x),np.array(a2F_x).shape)
        self.plot_lambda(Lambda_1)
        self.plot_lambda(self.lambda_2)
        lambda_2 = 2*np.trapz(self.lambda_2,dx=(Frequencies[9]-Frequencies[0])/10)
        #pass
        return Lambda_1, lambda_2
        #return Lambda_1

    def a2F_2(self,x):
        """
        This tries to calculate the Eliashberg function. It must be calculated in the full 'q'
        Inputs:
        x: The frequencies (but what frequencies, the plasmon ones?)
        """
        #print ("shape = ",np.shape(self.pars))
        center = self.pars[:,1]
        width = self.pars[:,2]
        #gauss_width = width/20
        #gauss_width = 0.0002 # what is the best width? It should be independent of the width to a point
        gauss_width = 0.01
        #print('gauss_width=',gauss_width)
        factor1 = 1/(2*self.Ne[np.where(self.energy==0.0)])
        #self.plot_lambda(center)
        #print('center=',center)
        center = np.absolute(center) #test to force the abs
        #self.plot_lambda(width)
        #print('width=',width)
        width = np.absolute(width)
        #self.plot_lambda(factor1)
        #print('factor1=',factor1)
        factor1 = np.absolute(factor1)
        suma = 0
        for i in range(len(center)):
            suma += (width[i]*center[i]) * self.gaussian(x,center[i],gauss_width) / (2*len(center)) #in qe-5.1.0_elph/elph_fc.f90 also has a weight related to Ne
                    #/len(center) # following the formula commented below # also choosing one element of the array
            # a2F(i) + gaussian * lambda(j,l) * w(j,l) * weightq(j) &
            #           / (2.0d0 * dble(total_qpoints_a2f))
            #print('width[',i,'],center[',i,']',width[i],center[i]) #OK
            print('gaussian(',x,',',center[i],')=', self.gaussian(x,center[i],gauss_width))
        #return factor1*suma
        #np.savetxt("a2F.txt",suma)
        #print('SUMA=',suma)
        return suma #I do not think a2F uses Ne factor

    def plot_lambda(self,x):
        plt.figure(figsize=(10,6))
        #plt.figure().set_figwidth(15)
        #plt.figure().set_figheight(2)
        plt.plot(x)
        plt.show()
        pass


def main(arg):
    #namefile = "xbv"
    #print (len(sys.argv))
    if (len( sys.argv ) > 1):
        print ("nombre de fichero:",arg[1])
        namefile = arg[1]
        print ("Creando el objeto plasmon")
        plasmon = Plasmon_analysis(arg,namefile)
        print ("Leyendo datos")

        if (len( sys.argv) == 3):
            index = arg[2]
            data, frequencies = plasmon.load_big_file(index, arg[1])
        else:
            data, frequencies, qx= plasmon.load_data()
#--------------------------diagonal-----------------------------
#        if (len(sys.argv[3]) == 4):
#            for i in range(51):
#                data[i,:] = plasmon.data[i,i,:]
#---------------------------------end---diagonal----------------
        print (np.shape(data),"=(51,5001)?")
        print ("dibuja")
        plasmon.plot_contour(data)
        print("calcula ajuste loretzian")
        #none = plasmon.locate_1Lorenztian(frequencies,data[30])
        #f=open('file_data_fittings.txt','a')
        #f.write("-amplitude--(+/-)error--center--(+/-)error--width--(+/-)error-\n")
        plasmon.fitting_Lorentz(frequencies,data)
        #f.close()
        if (len( sys.argv) == 3):
            np.savetxt("data_fitting_amplitudes_{}.txt".format(namefile),np.c_[plasmon.pars[:,0], plasmon.perr_lorentz[:,0]],header='#-----amplitude--(+/-)error---', footer='-------------')
            np.savetxt("data_fitting_center_{}.txt".format(namefile),np.c_[plasmon.pars[:,1], plasmon.perr_lorentz[:,1]],header='#-----center(e.V.)--(+/-)error---', footer='-------------')
            np.savetxt("data_fitting_width_{}.txt".format(namefile),np.c_[plasmon.pars[:,2], plasmon.perr_lorentz[:,2]],header='#-----width(e.V.)--(+/-)error---', footer='-------------')
        else:
            np.savetxt("data_fitting_amplitudes_{}.txt".format(namefile),np.c_[plasmon.pars[:,0], plasmon.perr_lorentz[:,0]],header='#-----amplitude--(+/-)error---for--qx={}'.format(qx), footer='-------------')
            np.savetxt("data_fitting_center_{}.txt".format(namefile),np.c_[plasmon.pars[:,1], plasmon.perr_lorentz[:,1]],header='#-----center(e.V.)--(+/-)error---for--qx={}'.format(qx), footer='-------------')
            np.savetxt("data_fitting_width_{}.txt".format(namefile),np.c_[plasmon.pars[:,2], plasmon.perr_lorentz[:,2]],header='#-----width(e.V.)--(+/-)error---for--qx={}'.format(qx), footer='-------------')

        print ("dibuja")
        plasmon.plot_contour(plasmon.Fitted_data)
        plasmon.plot_contour2(data,plasmon.Fitted_data)
        # superconductor = Eliashberg(plasmon.pars)
        # superconductor.read_Ne()
        ## lambda_1, lambda_2 = superconductor.Lambda(frequencies)
        # lambda_1 = superconductor.Lambda(frequencies)

        # print('Lambda_1=',np.sum(lambda_1)/len(lambda_1))
        # print('Lambda_2=',np.sum(superconductor.lambda_2))

        # np.savetxt('Lambda.txt', (lambda_1))
        # np.savetxt('Lambda_from_a2F.txt', np.array((frequencies[1:],superconductor.lambda_2)).T, header='frequencies,Lambda')
    else:
        print ("Arguments are namefile and the index of q_x as second argument if you want the BIG FILE")
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
