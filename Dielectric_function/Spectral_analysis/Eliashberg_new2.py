#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Eliashberg_new.py
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
# ---------------------------import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
from ase.units import create_units
import matplotlib.pyplot as plt
from scipy import integrate

import scipy as scipy
from scipy.signal import find_peaks
import os.path

import os
import glob
import xlrd

def excel_read(file_directory=".",header="1DP"):
    """
    Function that reads the data from the excel files in a directory.
    Header options are "1DP", "HPI", "HPII"
    """
    if (file_directory=="."):
        current_directory = os.getcwd()
    else:
        current_directory = file_directory
    # Find all .xlsx files in the current directory
    excel_files = glob.glob(os.path.join(current_directory, header+"*.xlsx"))
    data=[]
    # Iterate over each Excel file
    for file_name in excel_files:
        # Open the workbook
        workbook = xlrd.open_workbook(file_name)

        # Open the worksheet (assuming the first sheet is the one to be processed)
        worksheet = workbook.sheet_by_index(0)

        out=np.zeros((worksheet.nrows,worksheet.ncols))
        # Iterate the rows and columns
        for i in range(worksheet.nrows):
            row_values = []
            for j in range(worksheet.ncols):
                # Check if the cell is empty and fill with zero if true
                cell_value = worksheet.cell_value(i, j) if worksheet.cell_type(i, j) != xlrd.XL_CELL_EMPTY else  0
                # Add the cell value to the output array
                out[i][j] = cell_value

        data.append(out)
    #return qx,qy,out
    return data
def read_1_excel_file(file_directory=".",filename="1DP"):#filenames=('1DP','HPI','HPII')
    """
    This function only reads one excel file in qx,qy,w,e,r... format and outputs the data inside
    out=[[qx,qy,w,e,r,w2,e2,r2,w3,e3,r3],...]
    """
    if (file_directory=="."):
        current_directory = os.getcwd()
    else:
        current_directory = file_directory
    #excel_file = glob.glob(os.path.join(current_directory, filename+".xlsx"))
    # Open the workbook
    #workbook = xlrd.open_workbook(excel_file)
    workbook = xlrd.open_workbook(os.path.join(current_directory, filename+".xlsx"))
    # Open the worksheet (assuming the first sheet is the one to be processed)
    worksheet = workbook.sheet_by_index(0)

    out=np.zeros((worksheet.nrows,worksheet.ncols))
    # Iterate the rows and columns
    for i in range(worksheet.nrows):
        row_values = []
        for j in range(worksheet.ncols):
            # Check if the cell is empty and fill with zero if true
            cell_value = worksheet.cell_value(i, j) if worksheet.cell_type(i, j) != xlrd.XL_CELL_EMPTY else  0
            # Add the cell value to the output array
            out[i][j] = cell_value
    return out
    pass
def Excel_data_parser(data):
    qx=np.zeros(len(data))
    qy=np.zeros(len(data))
    Omega=np.zeros(len(data))
    Gamma=np.zeros(len(data))
    Ratio=np.zeros(len(data))
    for i in range(len(data)):
            qx[i]=data[i][0]
            qy[i]=data[i][1]
            Omega[i]=data[i][2]
            Gamma[i]=data[i][3]
            Ratio[i]=data[i][4]
    return qx,qy,Omega,Gamma,Ratio
    pass

    #Python program to reverse an array
def list_reverse(arr,size):

    #if only one element present, then return the array
    if(size==1):
        return arr

    #if only two elements present, then swap both the numbers.
    elif(size==2):
        arr[0],arr[1],=arr[1],arr[0]
        return arr

    #if more than two elements presents, then swap first and last numbers.
    else:
        i=0
        while(i<size//2):

    #swap present and preceding numbers at time and jump to second element after swap
            arr[i],arr[size-i-1]=arr[size-i-1],arr[i]

    #skip if present and preceding numbers indexes are same
            if((i!=i+1 and size-i-1 != size-i-2) and (i!=size-i-2 and size-i-1!=i+1)):
                arr[i+1],arr[size-i-2]=arr[size-i-2],arr[i+1]
            i+=2
        return arr
def plot_all(data,data2,pars,pars2,lambda_w_lista,w,a2F_lista,frequencies):
    """
    Plot all data
    """
    plt.rcParams['text.usetex'] = True
    fig = plt.figure(figsize=(10,6))
    ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2)
    im1 = ax1.imshow(data.T,
     #	vmin = 0.0 , vmax = 0.004,
     #	vmin = 0.0 , vmax = 0.175,
         vmin = 0.0 , vmax = 0.12,
         cmap=plt.colormaps['jet'], origin='lower',
         interpolation='gaussian', aspect='auto')
    ax1.set_title('Original data')
    ax1.set_xlabel('$q_x$ (a.u.)')
    ax1.set_ylabel('$\omega$ (eV)')
    ax1.set_xticks([0,51,101])
    ax1.set_yticks([0,5001])
    ax1.set_xticklabels(["($\pi$,0)","(0,0)","($\pi$,$\pi$)"])
    ax1.set_yticklabels(["0","1"])
    #
    ax2 = plt.subplot2grid((2, 3), (1, 0), colspan=2,  sharex=ax1, sharey=ax1)
    im2 = ax2.imshow(data2.T,
    #	vmin = 0.0 , vmax = 0.004,
    #	vmin = 0.0 , vmax = 0.3,
         vmin = 0.0 , vmax = 0.12,
         cmap=plt.colormaps['jet'], origin='lower',
         interpolation='gaussian', aspect='auto')
     #cax2.tick_params(direction='in', length=6, width=2, colors='k', right=True, labelright='on')
     #ax1[1].set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)
     #print (self.pars[:,1])
    x_max_range = len(pars[:,1])
    for xi in range(x_max_range):
        ax2.scatter(x = x_max_range-xi, y = pars[xi,1]*5001,c='k',marker='x',s=10)
        y_min =(pars[xi,1]-pars[xi,2])*5001
        y_max =(pars[xi,1]+pars[xi,2])*5001
        ax2.scatter(x = x_max_range-xi, y = y_min,c='k',marker='1',s=10)
        ax2.scatter(x = x_max_range-xi, y = y_max,c='k',marker='2',s=10)

        ax2.scatter(x = x_max_range+xi, y = pars2[xi,1]*5001,c='k',marker='x',s=10)
        y_min =(pars2[xi,1]-pars2[xi,2])*5001
        y_max =(pars2[xi,1]+pars2[xi,2])*5001
        ax2.scatter(x = x_max_range+xi, y = y_min,c='k',marker='1',s=10)
        ax2.scatter(x = x_max_range+xi, y = y_max,c='k',marker='2',s=10)
    ax2.axvline(x = len(pars[:,1]), linestyle = "--", color = "red")
    fig.colorbar(im1, ax=ax1)
    fig.colorbar(im2, ax=ax2)
#    cbar = fig.colorbar(ax1)
#    cbar2 = fig.colorbar(ax2)
    #fig.suptitle('$-Im(\epsilon^{-1}(q,\omega))$')
    #ax1.set_title('Original data')
    #ax1.set_ylabel('$\omega$ (eV)')
    ax2.set_ylabel('$\omega$ (eV)')
    #ax1.set_xlabel('$q_x$ (a.u.)')
    ax2.set_xlabel('$q_x$ (a.u.)')
    ax2.set_title('Lorentz fitting')
    #ax1.set_xticks([0,51,102])
    #ax1.set_yticks([0,5001])
    #ax1.set_xticklabels(["($\pi$,0)","(0,0)","($\pi$,$\pi$)"])
    #ax1.set_yticklabels(["0","1"])
    #
    ax4 = plt.subplot2grid((2, 3), (0, 2))#,  sharey=ax1)
    ax4.plot(w[1:],lambda_w_lista)
    ax4.set_title('$\lambda$ vs. $\omega$')
    ax4.set_ylabel('$\lambda$')
    ax4.set_xlabel('$\omega$ (eV)')
    ax5 = plt.subplot2grid((2, 3), (1, 2),  sharex=ax4)#, sharey=ax1)
    ax5.plot(frequencies,a2F_lista)
    ax5.set_title('$\\alpha^2F$ vs. $\omega$')
    ax5.set_xlabel('$\omega$ (eV)')
    ax5.set_ylabel('$\\alpha^2F$')
    #annotate_axes(fig)
    plt.tight_layout()
    plt.show()
    fig.savefig("Ajuste_total")

    #plt.show()
    return 0


class Plasmon_analysis(object):
    """Code for the analysis of Silkin's plasmon data."""

    def __init__(self, arg, namefile):
        super(Plasmon_analysis, self).__init__()
        self.arg = arg
        self.namefile = namefile
        self.pars = np.zeros((51,3))
        self.pars2 = np.zeros((51,3))
        self.perr_lorentz = np.zeros((51,3))
        self.perr_lorentz2 = np.zeros((51,3))

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
        self.data = np.resize(epsilon,(len(np. unique(q_x)),
                                        len(np. unique(q_y)),
                                        len(np. unique(w))))
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
        	vmin = 0.0 , vmax = 0.2,
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
        cax3 = ax1[1].scatter(x = range(len(self.pars[:,1])), y = (self.pars[:,1])*5001,c='k',marker='x',s=10)
        cax3 = ax1[1].scatter(x = range(len(self.pars[:,1])), y = (self.pars[:,1]+self.pars[:,2])*5001,c='k',marker='1',s=10)
        cax3 = ax1[1].scatter(x = range(len(self.pars[:,1])), y = (self.pars[:,1]-self.pars[:,2])*5001,c='k',marker='2',s=10)
        cax3 = ax1[1].axvline(x = len(self.pars[:,1]), linestyle = "--", color = "red")

        cbar = fig.colorbar(cax)
        cbar2 = fig.colorbar(cax2)
        if (len( sys.argv) == 3):
            plt.savefig("Ajuste_{}".format("original_vs_fitted"))
        else:
            plt.savefig("Ajuste_{}_{}".format(self.namefile,"original_vs_fitted"))
        plt.show()
        pass

    def plot_contour3(self,data,data2):
        plt.style.use('_mpl-gallery-nogrid')
        fig, ax1 = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)
        fig.set_size_inches(10, 5)
        fig.set_dpi(100)
        cax = ax1[0].imshow(data.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.175,
        	vmin = 0.0 , vmax = 0.12,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        #plt.tick_params(direction='in', length=6, width=2, colors='k', right=True, labelright='on')
        cax2 = ax1[1].imshow(data2.T,
        #	vmin = 0.0 , vmax = 0.004,
        #	vmin = 0.0 , vmax = 0.3,
            vmin = 0.0 , vmax = 0.12,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        x_max_range = len(self.pars[:,1])
        for xi in range(x_max_range):
            cax3 = ax1[1].scatter(x = x_max_range-xi, y = self.pars[xi,1]*5001,c='k',marker='x',s=10)
            y_min =(self.pars[xi,1]-self.pars[xi,2])*5001
            y_max =(self.pars[xi,1]+self.pars[xi,2])*5001
            cax3 = ax1[1].scatter(x = x_max_range-xi, y = y_min,c='k',marker='1',s=10)
            cax3 = ax1[1].scatter(x = x_max_range-xi, y = y_max,c='k',marker='2',s=10)

            cax3 = ax1[1].scatter(x = x_max_range+xi, y = self.pars2[xi,1]*5001,c='k',marker='x',s=10)
            y_min =(self.pars2[xi,1]-self.pars2[xi,2])*5001
            y_max =(self.pars2[xi,1]+self.pars2[xi,2])*5001
            cax3 = ax1[1].scatter(x = x_max_range+xi, y = y_min,c='k',marker='1',s=10)
            cax3 = ax1[1].scatter(x = x_max_range+xi, y = y_max,c='k',marker='2',s=10)
        cax3 = ax1[1].axvline(x = len(self.pars[:,1]), linestyle = "--", color = "red")

        cbar = fig.colorbar(cax)
        cbar2 = fig.colorbar(cax2)
        ax1[0].set_title('Original data')
        ax1[0].set_ylabel('$\omega$ (eV)')
        ax1[1].set_ylabel('$\omega$ (eV)')
        ax1[0].set_xlabel('$q_x$ (a.u.)')
        ax1[1].set_xlabel('$q_x$ (a.u.)')
        ax1[1].set_title('Lorentz fitting')
        ax1[0].set_xticks([0,51,102])
        ax1[0].set_yticks([0,5001])
        ax1[0].set_xticklabels(["($\pi$,0)","(0,0)","($\pi$,$\pi$)"])
        ax1[0].set_yticklabels(["0","1"])

        plt.tight_layout()
        if (len( sys.argv) == 3):
            plt.savefig("Ajuste_d_{}".format("original_vs_fitted"))
        else:
            plt.savefig("Ajuste_d_{}_{}".format(self.namefile,"original_vs_fitted"))

        plt.show()
        pass

    def locate_1Lorenztian(self, Frequency, Spec, index, big=False):
        def _1Lorentzian(x, amp1, cen1, wid1):
            return amp1*wid1**2/((x-cen1)**2+wid1**2)
        peaks, _ = find_peaks(Spec,width=5,rel_height=0.3)
        amp = 0.009
        cen = 0.009
        wid = 0.001
        popt_lorentz, pcov_lorentz = scipy.optimize.curve_fit(_1Lorentzian, Frequency, Spec, p0=[amp, cen, wid])

        perr_lorentz = np.sqrt(np.diag(pcov_lorentz))
        pars = popt_lorentz[:]
    #---------------------BIG---------------------------------------------------
        if big:
            if (Frequency[0] != 0):
                Frequency_big = np.append(Frequency,Frequency[:len(Frequency)//3]+Frequency[-1])
            else:
                Frequency_big = np.append(Frequency,Frequency[1:len(Frequency)//3]+Frequency[-1])
            lorentz_peak = _1Lorentzian(Frequency_big, *pars)
        else:
            lorentz_peak = _1Lorentzian(Frequency, *pars)
    #-------------------END_BIG-------------------------------------------------
        print ("-------------Peak ",index,"-------------")
        print ("amplitude = %0.2f (+/-) %0.2f" % (pars[0], perr_lorentz[0]))
        print ("center = %0.2f (+/-) %0.2f" % (pars[1], perr_lorentz[1]))
        print ("width = %0.2f (+/-) %0.2f" % (pars[2], perr_lorentz[2]))
        print ("area = %0.2f" % np.trapz(lorentz_peak))
        print ("--------------------------------")

        self.pars[index] = pars
        self.perr_lorentz[index] = perr_lorentz
        return lorentz_peak

    def fitting_Lorentz(self,frequencies,data,size,big=False):
        #self.Fitted_data = np.zeros((size,5001))
        self.Fitted_data = np.zeros((size,len(frequencies)))
        for i in range(size):
            self.Fitted_data[i] = self.locate_1Lorenztian(frequencies,data[i],i)
        #---------------------------new-----------------------------------------
        if big:
            if (frequencies[0] != 0):
                frequencies_big = np.append(frequencies,frequencies[:len(frequencies)//3]+frequencies[-1])
            else:
                frequencies_big = np.append(frequencies,frequencies[1:len(frequencies)//3]+frequencies[-1])
            self.Fitted_data = np.zeros((size,len(frequencies_big)))
            for i in range(size):
                self.Fitted_data[i] = self.locate_1Lorenztian(frequencies,data[i],i,big=True)
        #-----------------------------------------------------------------------
        pass

    def fitting_Lorentz2(self,frequencies,data,size,big=False):
        #self.Fitted_data2 = np.zeros((size,5001))
        self.Fitted_data2 = np.zeros((size,len(frequencies)))
        for i in range(size):
            self.Fitted_data2[i] = self.locate_1Lorenztian2(frequencies,data[i],i)
        #---------------------------new-----------------------------------------
        if big:
            if (frequencies[0] != 0):
                frequencies_big = np.append(frequencies,frequencies[:len(frequencies)//3]+frequencies[-1])
            else:
                frequencies_big = np.append(frequencies,frequencies[1:len(frequencies)//3]+frequencies[-1])
            self.Fitted_data2 = np.zeros((size,len(frequencies_big)))
            for i in range(size):
                self.Fitted_data2[i] = self.locate_1Lorenztian2(frequencies,data[i],i,big=True)
        #-----------------------------------------------------------------------
        pass

    def locate_1Lorenztian2(self, Frequency, Spec, index, big=False):
        def _1Lorentzian(x, amp1, cen1, wid1):
            return amp1*wid1**2/((x-cen1)**2+wid1**2)
        peaks, _ = find_peaks(Spec,width=5,rel_height=0.3)
        amp = 0.009
        cen = 0.009
        wid = 0.001
        popt_lorentz, pcov_lorentz = scipy.optimize.curve_fit(_1Lorentzian, Frequency, Spec, p0=[amp, cen, wid])

        perr_lorentz = np.sqrt(np.diag(pcov_lorentz))
        pars = popt_lorentz[:]
    #---------------------BIG---------------------------------------------------
        if big:
            if (Frequency[0] != 0):
                Frequency_big = np.append(Frequency,Frequency[:len(Frequency)//3]+Frequency[-1])
            else:
                Frequency_big = np.append(Frequency,Frequency[1:len(Frequency)//3]+Frequency[-1])
            lorentz_peak = _1Lorentzian(Frequency_big, *pars)
        else:
            lorentz_peak = _1Lorentzian(Frequency, *pars)
    #-------------------END_BIG-------------------------------------------------
        print ("-------------Peak ",index,"-------------")
        print ("amplitude = %0.2f (+/-) %0.2f" % (pars[0], perr_lorentz[0]))
        print ("center = %0.2f (+/-) %0.2f" % (pars[1], perr_lorentz[1]))
        print ("width = %0.2f (+/-) %0.2f" % (pars[2], perr_lorentz[2]))
        print ("area = %0.2f" % np.trapz(lorentz_peak))
        print ("--------------------------------")

        self.pars2[index] = pars
        self.perr_lorentz2[index] = perr_lorentz

        return lorentz_peak


#-------------------------------------------------------------
class Eliashberg2(object):
    """docstring for Eliashberg.
    This object calculates the Lambda from two methods:
    From the Eliashberg function amd fron the lambda at the q points.
    """

    def __init__(self, qx,qy,Omega,Gamma,Ratio):
        """
        Set parameters into object and transformation units.
        """
        units = create_units('2014')   #new way of defining units
        #self.pars = pars
        self.from_cm1_to_Hartree = units.invcm/units.Hartree #0.0000045563352812122295
        self.from_GHz_to_Hartree = self.from_cm1_to_Hartree /29.9792458 # from GHz to Hartree
        self.from_Ry_to_Hartree = units.Ry/units.Hartree #0.5 # from Ry to Hartree
        """
        cm^-1 to Hartree: Hartree = cm^-1 / 219474.6
        Hartree to cm^-1: cm^-1 = Hartree * 219474.6

        GHz to Hartree: Hartree = GHz / 6.5796839 × 10^9
        Hartree to GHz: GHz = Hartree * 6.5796839 × 10^9

        eV to Hartree: Hartree = eV*0,0367493
        """
        self.from_cm1_to_eV = units.invcm/units.eV #0.00012398425731484318
        self.from_GHz_to_eV = 0.000004135669661004052
        self.qx = qx
        self.qy = qy
        self.Omega = Omega
        self.Gamma =  Gamma
        self.Ratio = Ratio

        self.indice_zeros = 0
        self.N_qs = ((np.max(qx)*np.max(qy))+np.max([np.max(qx),np.max(qy)]))/(2*50)
        self.N = 50*50
        print("factor N_qs=",self.N_qs,"::",self.N,"::",len(Omega))

    def read_Ne(self,filename="out_DOS.dat"):
        """
        Read the number of particles under the Fermi level (???)
        Set the numbre of elements.
        """
        self.energy, self.Ne = np.loadtxt(filename,usecols=(0,1), unpack=True)
        numero_de_elementos = np.trapz(self.Ne[:(np.where(self.energy==0.0)[0][0])],dx=np.absolute(self.energy[0]-self.energy[-1])/len(self.energy)) # maybe this is better??
        numero_de_elementos2 = integrate.simpson(self.Ne[:(np.where(self.energy==0.0)[0][0])],dx = np.absolute(self.energy[0]-self.energy[-1])/len(self.energy))
        total_states = np.trapz(self.Ne,dx=np.absolute(self.energy[0]-self.energy[-1])/len(self.energy))
        print("Número de elementos:",numero_de_elementos,"|total=",total_states)
        print("Número de elementos2:",numero_de_elementos2,"|total=",total_states)
        print("N(eF)=",self.Ne[(np.where(self.energy==0.0)[0][0])]) #states/spin/eV/unit cell

        #------^--------
        print ("Density of states at Fermi level per cell=",self.Ne[np.where(self.energy==0.0)])
        print ("Number of elements (for the factor)=",numero_de_elementos)
        self.N_ef = self.Ne[(np.where(self.energy==0.0)[0][0])] #1.8855775
        pass

    def Lambda_q(self,gamma_q,omega_q,Nef):
        """
        Calculates the Lambda(q) functions
        ---input---
        gamma_q: widths of the lorentzian fitting of the plasmon
        omega_q: Frequencies of the plasmon, fitted to lorentzian
        Nef: Density of states at Fermi level per cell
        ---output---
        Lamb_q: Lambda(q)
        """
        Lamb_q=(1/(np.pi*Nef)) * (gamma_q/omega_q**2) #fix from omega to omega²
        return Lamb_q

    def Lambda_q_new(self,i):
        """
        Calculates the Lambda(q) functions
        ---input---
        i: index of the frequencies and widths of the lorentzian fitting of the plasmon
        ---output---
        Lamb_q: Lambda(q)
        """
        if (self.Omega[i]!=0):
            #print("OK_self.Gamma[",i,"]=",self.Gamma[i])
            Lamb_q=(1/(np.pi*self.N_ef)) * (self.Ratio[i])/(self.Omega[i]) #fix from omega to omega²
            # if Lamb_q<0:
            #     print("Lambda_q negativo??;",Lamb_q)
        else:
            print(self.Omega)
            print("self.Gamma[",i,"]=",self.Omega[i])
            #exit()
            Lamb_q=0
            self.indice_zeros += 1
            print("Indice_zeros=",self.indice_zeros)
        return Lamb_q

    def a2F(self,x):
        """
        Calculates the Eliashberg spectral functions
        ---input---
        x: coordiate to calculate the Eliashberg function.
        ---output---
        a2F = factor1*summa: The Eliashberg function at "x"
        """
        method=1
        center = self.Omega
        width = self.Gamma
        width = np.absolute(width)
        gauss_width = 0.004 #from 0.01 [0.1,0.05,0.01,0.005,0.001,0.0005,0.0001] 0.004 is the best option for the test
        summa = 0
        #---------method1---------------
        if (method==1):
            #factor1 = 1 / (2*(len(center)-self.indice_zeros)) #a2F(w)=1/2N Sum{Lambda*Omega*delta(w-Omega)}
            factor1 = 1 / (2*self.N) #a2F(w)=1/2N Sum{Lambda*Omega*delta(w-Omega)}
            for i in range(len(center)):
                summa += (width[i]*center[i]) * self.gaussian(x,center[i],gauss_width)
            return factor1*summa
        #--------method2------------------
        else:
            #factor2 = 1 / (2*np.pi*self.N_ef*(len(center-self.indice_zeros))) #a2F(w)=1/2piN(EF)N Sum{(Gamma/Omega)*delta(w-Omega)}
            factor2 = 1 / (2*np.pi*self.N_ef*self.N) #a2F(w)=1/2piN(EF)N Sum{(Gamma/Omega)*delta(w-Omega)}
            for i in range(len(center)):
                summa += (width/center) * self.gaussian(x,center[i],gauss_width)
            return factor2*summa

    def a2F_new(self,x,method = 0):
        """
        Calculates the Eliashberg spectral functions
        ---input---
        x: coordiate to calculate the Eliashberg function.
        method: Method to use in the calculation.
        ---output---
        a2F = factor1*summa: The Eliashberg function at "x"
        """

        center = self.Omega[:] #*put units correctly...
        width = self.Gamma[:] #*put units the same as center
        #width = np.absolute(width)
        gauss_width = 50000*self.from_cm1_to_eV#(units.invcm/units.Hartree) #0.00002 # test the units of this... should be aprox. 5 cm-1 (1, 5 or 10)
        if (method == 1):
            #---------method1-------vvvv---
            summa = 0
            for i in range(len(self.qx)):
                simmetry_factor =  2 #8
                if self.qx[i] ==  0:
                    simmetry_factor =  1 #4
                    if self.qy[i] ==  0:
                        simmetry_factor =  0 #1
                elif self.qy[i] ==  0:
                    simmetry_factor =  1 #4
                elif self.qx[i] == self.qy[i]:
                    simmetry_factor =  1 #4
                summa += self.Lambda_q_new(i)*self.Omega[i] * self.gaussian(x,self.Omega[i],gauss_width)
            #return summa/(2*(len(self.Omega)-self.indice_zeros))
            return summa/(2*self.N)
        else:
            #---------method2------vvvv-
            summa = 0

            for i in range(len(self.qx)):
                simmetry_factor =  2 #8
                if self.qx[i] ==  0:
                    simmetry_factor =  1 #4
                    if self.qy[i] ==  0:
                        simmetry_factor = 0 #1
                elif self.qy[i] ==  0:
                    simmetry_factor =  1 #4
                elif self.qx[i] == self.qy[i]:
                    simmetry_factor =  1 #4

                summa += (self.Ratio[i]) * self.gaussian(x,self.Omega[i],gauss_width) * simmetry_factor
            #return summa/(2*np.pi*self.N_ef*(len(self.Omega)-self.indice_zeros))
            #return summa/(2*np.pi*self.N_ef*len(self.Omega))
            #return summa/(2*np.pi*self.N_ef*self.N_qs)
            return summa/(2*np.pi*self.N_ef*self.N)

    def gaussian(self,x, center,gauss_width):
        """
        Gaussian distribution
        ---input---
        x: x position of the gaussian graph, in this case the frequency
        center: center of the gaussian distribution, in this case the frequency of the plasmon Lorenztian approximation
        gauss_width: width of the gaussian distribution, in this case the witdh of the plasmon lorentzian aproximation
        ---output---
        g: gaussian value in 'x'
        """
        mu = center #self.pars[:1] #center
        sigma = gauss_width #self.pars[:,2] #width
        #d = 1/(2*np.sqrt(2*np.pi)*sigma)
        d = 1/(sigma*np.sqrt(2*np.pi))
        g = d*np.exp(-(x-mu)**2 / ( 2.0 * sigma**2 )  )
        #exp(-(w_aux-w(j,l))**2.0d0/(2.0d0*broad**2.0d0))/ (broad*sqrt(twopi)) #copy from qe-5.1.0_elph
        return g

    def test_gaussian(self,w,w_q,N_q):
        """
        Test for the gaussian, the result must be the DOS.
        ---input---
        w: frequencies to test
        w_q: frequencies of the plasmon grid
        N_q: number of terms in w_q
        ---output---
        test_dos:
        """
        width = 5#*self.from_cm1_to_eV #self.from_cm1_to_Hartree*2 #1,5,10 cm-1
        test_dos = 1/N_q # is this equal to len(W_q)?
        for i in range(len(w_q)): #q=6x6x6 (???) for the test data. q=50x50=250 for the plasmon
            test_dos += self.gaussian(w,w_q[i],width) #This must be the DOS, integrate this to check the gaussian.
        print("Test of the gaussian: this must be the integral of the DOS=",test_dos)
        return test_dos


    def plot_lambda(self,x,w):
        plt.figure(figsize=(10,6))
        #plt.figure().set_figwidth(15)
        #plt.figure().set_figheight(2)
        plt.plot(x,w)
        plt.show()
        pass

    def Lambda(self,Frequencies):
        """
        Calculates the Lambda by two methods, notice that it must calculate the integral
        in a range that takes the Lorenztian obtained by the Plasmon_analysis object.
        """
        center = self.Omega
        width = self.Gamma
        #print('len(Center)-indice_zeros=',len(center)-self.indice_zeros)
        #Nef = self.Ne[np.where(self.energy==0.0)[0][0]]  #test
        Nef = self.N_ef
        center = np.absolute(center) #test to force the abs
        width = np.absolute(width)
        summa1 = 0
        for i in range(len(center)):
            summa1 += self.Lambda_q(width[i],center[i],Nef)
        #Lambda_1=summa1/(len(center)-self.indice_zeros) #* 33 #misterious factor... joking, this is the number of nodes in the example.
        Lambda_1=summa1/self.N
        self.lambda_2=[]
#        Frequncies = Frequencies+1
#        Frequencies = np.append(Frequencies,Frequncies, axis=0)
        if (Frequencies[0] != 0):
             Frequencies = np.append(Frequencies,Frequencies[:len(Frequencies)//3]+Frequencies[-1])
        else:
             Frequencies = np.append(Frequencies,Frequencies[1:len(Frequencies)//3]+Frequencies[-1])
#        for w in Frequencies[0:int(2*len(Frequencies)/3)]:
        for w in Frequencies[0:int(len(Frequencies))]: #edit for test
             print('Frequencie(',w,')              ', end="\r", flush=True)
             if (w == 0):
                 a2F_x = []
             else:
                 a2F_x.append(self.a2F(w)/w)
#---test--v--multiprocessing*****
        w = Frequencies[Frequencies != 0]
        # w = np.append(w,w[:len(w)//3]+w[-1]) # to expand the frequencies to cover the widths (and some extra) for the calculations
        with mp.Pool() as pool:
           res = pool.map(self.a2F,w)
        # a2F_x = res / w #or np.divide(res, w)
        a2F_x = np.divide(res, w)
#---test--^--multiprocessing*****
        #self.lambda_2 = 2*np.trapz(a2F_x,dx=(Frequencies[-1]-Frequencies[0])/len(Frequencies)) <-- this makes lambda1 20 times bigger than lambda2 (???)
        self.lambda_2 = 2*np.trapz(a2F_x) #test <-- this lambda2 is 250 times lambda1 (that corresponds to q_x*q_y ???)
        self.plot_lambda(a2F_x,w)
        return Lambda_1

    def Lambda_new(self,Frequencies):
        """
        Calculates the Lambda by two methods, notice that it must calculate the integral
        in a range that takes the Lorenztian obtained by the Plasmon_analysis object.
        """
        center = self.Omega
        width = self.Gamma
        #print('len(Center)-indice_zeros=',len(center)-self.indice_zeros)
        #method 1 ----------------------------

        summa1 = 0
        self.lambda_q_lista = np.array([])
        for i in range(len(center)):  #summa in q (6x6x6 or q_x*q_y)
            simmetry_factor =  2 #8
            if self.qx[i] ==  0:
                simmetry_factor =  1 #4
                if self.qy[i] ==  0:
                    simmetry_factor =  0 #1 #test taking out the qx=qy=0
            elif self.qy[i] ==  0:
                simmetry_factor =  1 #4
            elif self.qx[i] == self.qy[i]:
                simmetry_factor =  1 #4

            summa1 += self.Lambda_q_new(i)*simmetry_factor
            self.lambda_q_lista = np.append(self.lambda_q_lista, self.Lambda_q_new(i))
        #Lambda_1=summa1/(len(center)-self.indice_zeros)
        Lambda_1=summa1/self.N
        print("Test Lambda with lambda_q=",summa1/(len(center)),":test lambda_q[N]=",summa1/self.N,":lambda_q[N_qs]=",summa1/self.N_qs)
        #method 2 -------------------------------
        self.lambda_2=[]
        w = Frequencies[Frequencies != 0]
        self.w = w
        #w = np.linspace(0.0001,0.9999,20000) #test
        # with mp.Pool() as pool:
        #     res = pool.map(self.a2F_new,w)
        res = self.a2F_new(w)
        a2F_x = np.divide(res, w)
        #self.lambda_2_test2 = 2*integrate.simpson(np.divide(res,w)*self.from_cm1_to_Hartree /29.9792458,w)
        #self.lambda_2 = 2*integrate.simpson(a2F_x,w)
        self.lambda_2 = 2*np.trapz(a2F_x, w)
        print("plot a2F(w)")
        self.plot_lambda(res,w)
        print("plot Lambda(w)")
        #self.lambda_w_lista = []
        #-------removin negatives w's---
        mask = w >=  0
        self.w_0 = w[mask]
        self.lambda_w_lista = np.zeros(len(self.w_0))
        #print("a2F(",i,")/w=",np.divide(self.a2F_new(self.w_0), self.w_0)) #test
        diff_w_0 = np.diff(self.w_0)
        print("w_0 always growing?=",np.all(diff_w_0 > 0)) # this test is OK
        for i in range(1,len(self.w_0)): #Frequencies[Frequencies != 0]:
            #print("a2F(",i,")=",self.a2F_new(w[i])) #test = all this is > 0
            #print("a2F(",i,")/w=",np.divide(self.a2F_new(self.w_0[i]), self.w_0[i])) #test = all this is > 0
            w_1 = self.w_0[:i]
            #func_w = np.divide(self.a2F_new(w_1), w_1)
            func_w = np.absolute(self.a2F_new(w_1)/w_1)
            # mask = w_1 <= 0
            # if len(w_1[mask])>0:
            #     print("w_1 error")
            # mask = func_w <=0
            # if len(func_w[mask])>0:
            #     print("func_w error")
            # if len(func_w)!=len(w_1):
            #     print("error:",len(func_w),"::",len(w_1))
            #partial_value=2*integrate.simpson(func_w,w_1)
            partial_value = 2* np.trapz(func_w,w_1)
            #print("Lambda[",i,"]=",partial_value)
            #self.lambda_w_lista.append(partial_value) #test to see the evolution of Lambda
            self.lambda_w_lista[i] = partial_value #test to see the evolution of Lambda
        #self.plot_lambda(self.lambda_w_lista,self.w_0[1:])
        diff_lambda_w = np.diff(self.lambda_w_lista)
        print("lambda(w) always growing?=",np.all(diff_lambda_w > 0)) # this test fails!!!
        self.plot_lambda(self.lambda_w_lista,self.w_0)
        #---------end of removing negatives test-----
        # for i in range(1,len(w)): #Frequencies[Frequencies != 0]:
        #     #print("a2F(",i,")=",self.a2F_new(w[i])) #test = all this is > 0
        #     if w[i]>0:
        #         print("a2F(",i,")/w=",np.divide(self.a2F_new(w[i]), w[i])) #test
        #     w_1 = w[:i]
        #     partial_value=2*integrate.simpson(np.divide(self.a2F_new(w_1), w_1),w_1)
        #     #print("Lambda[",i,"]=",partial_value)
        #     self.lambda_w_lista.append(partial_value) #test to see the evolution of Lambda
        # self.plot_lambda(self.lambda_w_lista,w[1:])
        return Lambda_1


    def T_c(self,mu_par,lambda_t):
        """
        Allen-Dynes formula
        Tc=W_log/1.2*exp(-1.04*(1+self.lambda_2)/(self.lambda_2-mu_*(1+0.62*self.lambda_2)))
        """
        AllenDynes = True
        # out2=-1.04*(1+self.lambda_2)/(self.lambda_2-mu_par*(1+0.62*self.lambda_2))
        out2=-1.04*(1+lambda_t)/(lambda_t-mu_par*(1+0.62*lambda_t))
        if (AllenDynes):
            #out_AD=f_1(mu_par)*f_2(mu_par)*self.w_log()/1.2 * np.exp(out2)
            print ("f1,f2=",self.f_1(mu_par,lambda_t),self.f_2(mu_par))
            out=self.f_1(mu_par,lambda_t)*self.f_2(mu_par)*self.w_log(lambda_t)/1.2 * np.exp(out2)
        else:
            out=self.w_log(lambda_t)/1.2 * np.exp(out2)

        return out

    def w_log(self,lambda_t):
        """
        w_log calculated from Eliashberg
        """
        eV_to_K=11604
        #eV_to_K = 11604.5250061657
        print("w=",self.w)
        mask = self.w >=  0
        w = self.w[mask]#*eV_to_K
        # w_log = np.exp(2.0/self.lambda_2*
        w_log = np.exp(2.0/lambda_t*
            #integrate.simpson(
            #(np.divide(self.a2F_new(self.w), self.w)*np.log(self.w)),self.w))
            np.trapz(
            (np.divide(self.a2F_new(w), w)*np.log(w)),w))
        return w_log

    def w_2(self,lambda_t):
        eV_to_K=11604
        w = self.w#*eV_to_K
        # print ("**********************************************************")
        # print ("w=",w)
        # print ("a2F_new(w)=",np.nonzero(self.a2F_new(w)))
        # self.plot_lambda(self.a2F_new(w))
        # print ("self.a2F_new(w)*w=",self.a2F_new(w)*w)
        # self.plot_lambda(self.a2F_new(w)*w)
        # print ("**********************************************************")

        # w_2 = np.sqrt(2.0/self.lambda_2*
        w_2 = np.sqrt(2.0/lambda_t*
            #integrate.simpson(
            np.trapz(
            (self.a2F_new(w)*w),w))
            #np.prod([self.a2F_new(w),w]),w))
        return w_2

    def f_1(self,mu_par,lambda_t):
        LAMBDA_temp = 2.46*(1+3.8*mu_par)
        # return np.power(1+np.power(self.lambda_2/LAMBDA_temp,3/2),1/3)
        return np.power(1+np.power(lambda_t/LAMBDA_temp,3/2),1/3)

    def f_2(self,mu_par,lambda_t):
        LAMBDA_temp =1.82*(1+6.3*mu_par)*(self.w_2(lambda_t)/self.w_log(lambda_t))
        print("mu_par=",mu_par)#OK
        print ("LAMBDA_temp=",LAMBDA_temp) #ERROR
        print ("w_2(lambda)=",self.w_2(lambda_t)) #ERROR
        print ("w_log(lambda)=",self.w_log(lambda_t))
        # print ("return f2=",1 + (( self.w_2()/self.w_log() - 1) * self.lambda_2**2)/(
        # (self.lambda_2**2) + (LAMBDA_temp**2)))
        # return 1 + (( self.w_2()/self.w_log() - 1) * self.lambda_2**2)/(
        # (self.lambda_2**2) + (LAMBDA_temp**2))
        print ("return f2=",1 + (( self.w_2(lambda_t)/self.w_log(lambda_t) - 1) * lambda_t**2)/(
        (lambda_t**2) + (LAMBDA_temp**2)))
        return 1 + (( self.w_2(lambda_t)/self.w_log(lambda_t) - 1) * lambda_t**2)/(
        (lambda_t**2) + (LAMBDA_temp**2))

    def plot_contour_d(self,data,mask_value=50, diagonal = True):
        print ("Leyendo datos")
        plasmon = Plasmon_analysis(4,"A7_EPS.dat")
        if diagonal:
            mask = self.qx == self.qy
            data, frequencies = plasmon.load_big_file(4, diagonal=True)
        else:
            mask = self.qx == mask_value
            data, frequencies = plasmon.load_big_file(4, diagonal=False)

        plt.style.use('_mpl-gallery-nogrid')
        fig, ax1 = plt.subplots(1,1)
        fig.set_size_inches(10, 5)
        fig.set_dpi(100)
        cax = ax1.imshow(data.T,
        #	vmin = 0.0 , vmax = 0.004,
        	vmin = 0.0 , vmax = 0.2,
        	cmap=plt.colormaps['jet'], origin='lower',
        	interpolation='gaussian', aspect='auto')
        ax1.set_ylabel(r'Frequency (cm$^{-1}$)', fontsize=12)


        print (mask)
        cax2 = ax1.scatter(x = self.qy[mask]/50, y = self.Omega[mask]*5,c='k',marker='x',s=10)
        cax2 = ax1.scatter(x = self.qy[mask]/50, y = (self.Omega[mask]+self.Gamma[mask]/2)*5,c='k',marker='1',s=10)
        cax2 = ax1.scatter(x = self.qy[mask]/50, y = (self.Omega[mask]-self.Gamma[mask]/2)*5,c='k',marker='2',s=10)
        cax3 = ax1.errorbar(x = self.qy[mask]/50, y = self.Omega[mask]*5, yerr = self.Gamma[mask]*5/2, fmt = 'o') #¿y*5001? para la escala
        #cbar = fig.colorbar(cax)
        print("qx=",self.qx[1:10])
        plt.show()

        pass
#-------------------------------------------------------------------------------
def main(arg):

    if (len( sys.argv ) > 1):
        print ("nombre de fichero:",arg[1])
        namefile = arg[1]
        print ("Creando el objeto plasmon")
        plasmon = Plasmon_analysis(arg,namefile)
        print ("Leyendo datos")

        if (len( sys.argv) == 3):
            index = arg[2]

            data, frequencies = plasmon.load_big_file(index, arg[1],diagonal=False)
            data_d, frequencies_d = plasmon.load_big_file(index, arg[1], diagonal=True)
            #data.append(data_d)
        #    data = np.vstack((np.flip(data, axis=1), data_d))
            #frequencies.append(frequencies_d)
        #    frequencies = np.vstack((np.flip(frequencies), frequencies_d))
        else:
            data, frequencies, qx= plasmon.load_data()
#--------------------------diagonal-----------------------------
#        if (len(sys.argv[3]) == 4):
#            for i in range(51):
#                data[i,:] = plasmon.data[i,i,:]
#---------------------------------end---diagonal----------------
#-------------------------fix to remove extra polaron data--------------
        if (len( sys.argv) == 3):
            for i in range(5):
                for j in range(2000,5001,1):
                    if (data[i][j]!=0.0):
                        print("data[{}][{}]=".format(i,j),data[i][j])
            # if i < 5:
            #     if j > 2000:
                        data[i][j] = 0
                    if (data_d[i][j]!=0.0):
                        print("data_d[{}][{}]=".format(i,j),data_d[i][j])
            # if i < 5:
            #     if j > 2000:
                        data_d[i][j] = 0
        else:
            for i in range(5):
                for j in range(2000,5001,1):
                    if (data[i][j]!=0.0):
                        print("data[{}][{}]=".format(i,j),data[i][j])
            # if i < 5:
            #     if j > 2000:
                        data[i][j] = 0
#------------------------end-removing data------------------------------
        print ("data[49][2110]=",data[49][2110]) # aqui hay datos
        print ("data[50][2110]=",data[50][2110]) # aqui no hay datos???
        print("***********************")
        print (np.shape(data),"=(51,5001)?")
        print ("Frequencies length=",len(frequencies))
        print ("dibuja")
        plasmon.plot_contour(data)
        #plasmon.plot_contour2(data,plasmon.Fitted_data)
    else:
        print ("Arguments are namefile and the index of q_x as second argument if you want the BIG FILE")
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
