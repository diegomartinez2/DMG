#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  TODO_multitherm.py (Version 0.0.2)
#
#  Copyright 2017 Diego Martinez Gutierrez
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

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import heapq
class data(object):
    def __init__(self):
        #Recuerda cambiar el 'timestep =' de forma adecuada.
        self.timestep = 0.001    #paso de tiempo del calculo
        #Area = 1.54*19.9*1e-20  #grosor X ancho en metros
        self.Area = 1.54*24.24*1e-20  #grosor X ancho en metros
        self.Cadenas = 3    #numero de cadenas moleculares
        #    recuerda cambiar el 'indice =' de forma adecuada.
        self.indice = range(25,122+1)+range(124,247+1)+range(249,345+1) #indice de los atomos que se mueven y que no son del termostato
        self.factorstep = 0            #Recuerda cambiar esto, el factor de desplazamiento hacia atras del punto final (para mirar la evolucion)
        self.lenghtinterval = 170000        #y esto es el numero de pasos de calculo desde el final
        self.num_pasos_tot = int(raw_input("Number of steps in the MD simulation?: "))
#        self.nombre_fichero_in = raw_input("Name of the file with the 'xyz' data: ")
#        self.nombre_fichero_in_2 = raw_input("Name of the file with the Molecular dynamics ('md.out') data: ") #This usually does not change
#        self.nombre_fichero_out = raw_input("Name of the output file: ")
#        if nombre_fichero_in == nombre_fichero_out:
#            raise ValueError("File input and output are the same!." )
#            return 1

class dTdx (object):

    def __init__(self):
        self.T_atoms = []
        self.Coordinates_X = []
        self.Coordinates_Y = []
        self.Coordinates_Z = []
        self.num_atoms = 0
        self.dTdx = 0
        self.DeltaT = 0

#        self.nombre_fichero_in = raw_input("Name of the file with the data: ")
#        self.nombre_fichero_out = raw_input("Name of the output file: ")
#        if self.nombre_fichero_in == self.nombre_fichero_out:
#            raise ValueError("File input and output are the same!." )
#            return 1

    def read_xyz(self,filename_in,filename_out,datos):
        """Lee las coordenadas de un archivo XYZ. Incluso las velocidades de acuerdo con el fichero geo_out.gen.
        Si el numero de coordenadas no coincide con lo establecido en el fichero da un aviso de error.
        Diego Martinez Gutierrez"""

        self.xyz = open(filename_in,"r")
        self.xyz_out = open(filename_out,"w")
        printProgress(0,datos.num_pasos_tot, prefix = 'Loading data:', suffix = ' ', barLength = 50)
        for tiempo in range(datos.num_pasos_tot):
            if (float(tiempo)/100.0) == float(tiempo/100): #progress bar no es necesario con pocos datos
                printProgress(tiempo,datos.num_pasos_tot, prefix = 'Loading data:', suffix = ' ', barLength = 50)
            self.num_atoms = int(self.xyz.readline().split()[0])
            title = self.xyz.readline()
            MD_step = title.split()[2]

            texto = str(self.num_atoms),"\n" #Number of atoms
            self.xyz_out.write(" ".join(texto))                    #Numero de atomos
            self.xyz_out.write("XYZ with local Temperatures\n")    #Just text
            self.xyz_out.write(title)
            for step in range(self.num_atoms):
                atom,x,y,z,vx,vy,vz = self.xyz.readline().split()
                T_tot,T=self.Temperature_local(atom,vx,vy,vz)
                if atom == 'C':
                    self.T_atoms.append(float(T))
                    if tiempo == 0:
                        self.Coordinates_X.append(float(x))
                        self.Coordinates_Y.append(float(y))
                        self.Coordinates_Z.append(float(z))
                texto = '{:5}{:15}{:15}{:15}{:15}{:15}{:15}{:15}{:15}{}'.format(atom,str(x),str(y),str(z),str(vx),str(vy),str(vz),str(T_tot),str(T),"\n")
                self.xyz_out.write("".join(texto))

        self.xyz.close()
        self.xyz_out.close()

        return 0

    def m_atomo(self,atom):

        lista_A= [['C',12.01],['H',1.008],['O',16.0]] #"""relative atomic masses (x1,660538921e−27 kg)"""
        for r in range(len(lista_A)):
            if lista_A[r][0]==atom:
                retorno = ((lista_A[r][1])*1.660538921e-27) #"""atomic mass in Kg"""
                return  retorno
        print ("atom not identified")
        return 1

    def Temperature_local(self,atom,vx,vy,vz):
        Boltzmann = 0.68749e+27
        Grados_de_libertad = 3
        T_tot = 0

        T = ((self.m_atomo(atom)*(float(vx)**2+float(vy)**2+float(vz)**2))/(Grados_de_libertad))*Boltzmann
        T_tot += T

        return T_tot, T


    def printGraphic(self,datos):
        x = []
        y = []
        z = []
        t = []
        t2 = []
        printProgress(0,datos.num_pasos_tot, prefix = 'Plotting data:', suffix = ' ', barLength = 50)
        for i_i in range(datos.num_pasos_tot):
            if (float(i_i)/1000.0) == float(i_i/1000): #progress bar no es necesario con pocos datos
                printProgress(i_i,datos.num_pasos_tot, prefix = 'Plotting data:', suffix = ' ', barLength = 50)
            for i_j in range(len(self.Coordinates_X)):
                if i_i == 0:
                    x.append(self.Coordinates_X[i_j])
                    y.append(self.Coordinates_Y[i_j])
                    z.append(self.Coordinates_Z[i_j])
                    t.append(self.T_atoms[i_j])
                if i_i != 0:
                    t[i_j] += self.T_atoms[(i_i*len(self.Coordinates_X))+i_j]
                    if i_i == (datos.num_pasos_tot)-datos.lenghtinterval*(1+datos.factorstep):
                        t2.append(self.T_atoms[(i_i*len(self.Coordinates_X))+i_j])
                    if ((datos.num_pasos_tot)-(datos.lenghtinterval*(1+datos.factorstep)) < i_i < (datos.num_pasos_tot)-(datos.lenghtinterval*datos.factorstep)):
                        t2[i_j] += self.T_atoms[(i_i*len(self.Coordinates_X))+i_j] # comprobar que todo esto esta bien ??? creo que es OK
    #--------------------X-T-------------------------
    #We renormalize the temperature and write usefull data to 'workfile'
        f = open('workfile','w')
        f2 = open('workfile2','w')
        for i_i in range(len(t)):
            t[i_i]/=float(datos.num_pasos_tot)
            t2[i_i]/=float(datos.lenghtinterval)
            if y[i_i] >= 23.0 and y[i_i] <= 27.0:
                f.write('{} {} {}\n'.format(x[i_i],t[i_i],t2[i_i]))
            if i_i <=36 :
                f2.write('{} {} {}\n'.format(x[i_i],t[i_i],t2[i_i]))
        f.close()
        f2.close()
    #-----------------------------------------------
        xmin = min(x)
        xmax = max(x)
        ymin = min(y)
        ymax = max(y)
    #    ymin = min(z)
    #    ymax = max(z)
    #    y = z
    #---------------------------------
        plt.figure(1)
    #+++++++++++++++++++++++++++++++++
    #    plt.subplot(211)
        plt.subplot(111)
        # define grid.
        xi = np.linspace(xmin,xmax, 200)
        yi = np.linspace(ymin,ymax, 200)
        # grid the data.
        for ii in range(len(t2)):
            try:
                del x[t2.index(0)]
                del y[t2.index(0)]
                del z[t2.index(0)]
                del t2[t2.index(0)]
            except:
                pass

    #    zi = griddata(x, y, t, xi, yi)#, interp='linear')
        zi = griddata(x, y, t2, xi, yi)#, interp='linear')
        # contour the gridded data, plotting dots at the nonuniform data points.
        CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
    #    CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow,vmax=abs(zi).max(), vmin=abs(zi).min())
    #    CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow,vmax=max(t), vmin=min(t))
        CS = plt.contourf(xi, yi, zi, 30, cmap=plt.cm.rainbow,vmax=max(t2), vmin=heapq.nsmallest(5,t2)[-1])
        plt.colorbar()  # draw colorbar
        # plot data points.
        plt.scatter(x, y, marker='o', c='b', s=5, zorder=10)
        plt.xlim(xmin,xmax)
        plt.ylim(ymin,ymax)
        plt.title('griddata of local Temperatures')
        plt.xlabel('$\AA$')
        plt.ylabel('$\AA$')
    #++++++++++++++++++++++++++++++++++
    #    plt.subplot(212) #the other plot
        #define the second grid
    #    xi2 = np.linspace(x2min,x2max, 200)
    #    yi2 = np.linspace(y2min,y2max, 200)
        # grid the data.
    #    zi2 = griddata(x2, y2, t22, xi2, yi2)#, interp='linear')
        # contour the gridded data, plotting dots at the nonuniform data points.
    #    CS = plt.contour(xi2, yi2, zi2, 15, linewidths=0.5, colors='k')
    #    CS = plt.contourf(xi2, yi2, zi2, 15, cmap=plt.cm.rainbow,vmax=abs(zi2).max(), vmin=abs(zi2).min())
    #    plt.colorbar()  # draw colorbar
        # plot data points.
    #    plt.scatter(x2, y2, marker='o', c='b', s=5, zorder=10)
    #    plt.xlim(x2min,x2max)
    #    plt.ylim(y2min,y2max)
    #    plt.title('griddata of local Temperatures of the device')
    #+++++++++++++++++++++++++++++++++++
        plt.show()
    #-----------full:
        x_1=x[:]
        t_1=t2[:]
        for ii in range(len(x_1)):
            if ii==0:
                x_2=[]
                t_2=[]
            if t_1[ii] > 0 and (0 < x_1[ii] < 200):
                x_2.append(x_1[ii])
                t_2.append(t_1[ii])
        Lin_regres(x_2,t_2,2,0)
    #----------and the parts:
        x_1=x[24:107]
        t_1=t2[24:107]
        for ii in range(len(x_1)):
            if ii==0:
                x_2=[]
                t_2=[]
            if t_1[ii] > 0 and (0 < x_1[ii] < 200):
                x_2.append(x_1[ii])
                t_2.append(t_1[ii])
        m_1,b_1=Lin_regres(x_2,t_2,0,1)
        x_1=x[-107:-24]
        t_1=t2[-107:-24]
        for ii in range(len(x_1)):
            if ii==0:
                x_2=[]
                t_2=[]
            if t_1[ii] > 0 and (0 < x_1[ii] < 200):
                x_2.append(x_1[ii])
                t_2.append(t_1[ii])
        m_2,b_2=Lin_regres(x_2,t_2,0,1)
        x_1=x[121:-121]
        t_1=t2[121:-121]
        for ii in range(len(x_1)):
            if ii==0:
                x_2=[]
                t_2=[]
            if t_1[ii] > 0 and (0 < x_1[ii] < 200):
                x_2.append(x_1[ii])
                t_2.append(t_1[ii])
        Lin_regres(x_2,t_2,1,4)
        self.dTdx = (abs(m_1)+abs(m_2))/2.0
        self.DeltaT = abs((m_1*110+b_1)-(m_2*110+b_2))
        print('media dT/dx= {} Kelvin/Armstrong ;error= {}'.format(self.dTdx,(abs(m_1)-abs(m_2))/2.0))
        print('Delta T= {} Kelvin'.format(self.DeltaT))
        print('verificando: {} en 102:: {} en 119'.format((m_1*102+b_1),(m_2*119+b_2)))
        plt.show()

        return 0

# Print iterations progress
def printProgress (iteration, total, prefix = '', suffix = '', decimals = 1, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    formatStr = "{0:." + str(decimals) + "f}"
    percent = formatStr.format(100 * (iteration / float(total)))
    filledLength = int(round(barLength * iteration / float(total)))
    bar = '+' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percent, '%', suffix)),
    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()



def Lin_regres(x,y,i,j):
    m,b = np.polyfit(x,y,1)
    fit = np.polyfit(x,y,j)
    fit_fn = np.poly1d(fit)
    if i == 0:
        plt.plot(x,y, 'yo', x, fit_fn(x), '--k')
        print(m,b,fit)
    if i == 1:
        plt.plot(x,y, 'ro', x, fit_fn(x), '--k')
    if i == 2:
        plt.plot(x,y, 'ro')
    plt.xlabel('$\AA$')
    plt.ylabel('$K$')
    return m,b


class dQdt(object):

    def __init__(self):
        self.dQdt = 0
        self.Heat1 = 0.0
        self.Heat2 = 0.0
        self.Heat3 = 0.0

    def read_md(self,filename_in,Heat1,Heat2,Heat3,datos):
        """Reads the heat exchange from md.out file and add the Heat exchanged"""
        self.Heat_1=[]
        self.Heat_3=[]
        with open(filename_in,"r") as md_file:
            while True:
                line_md = md_file.readline().split()
                if not line_md: break
                if line_md[0]=="Heat":
                    if line_md[2]=="(1):":
                        self.Heat1 += float(line_md[3])
                        self.Heat_1.append(self.Heat1)
                    if line_md[2]=="(2):":
                        self.Heat2 += float(line_md[3])
                    if line_md[2]=="(3):":
                        self.Heat3 += float(line_md[3])
                        self.Heat_3.append(self.Heat3)
        plt.plot(self.Heat_1[-datos.lenghtinterval::1],self.Heat_3[-datos.lenghtinterval::1])
        plt.show()
        x2=[x*datos.timestep  for x in range(len(self.Heat_1))]
        plt.plot(x2,self.Heat_1,x2,self.Heat_3)
        plt.xlabel('picoseconds')
        plt.ylabel('Hartree')
        plt.show()
        self.plot_Heats(x2,datos)
        m_1,b_1=Lin_regres(x2[-datos.lenghtinterval::1],self.Heat_1[-datos.lenghtinterval::1],0,1)
        m_2,b_2=Lin_regres(x2[-datos.lenghtinterval::1],self.Heat_3[-datos.lenghtinterval::1],0,1)
        self.dQdt = ((m_1+m_2)/2.0)*4.359745e-6
        print('dQ/dt={} Hartree/ps;={} W ;error={} W'.format((m_1+m_2)/2.0,self.dQdt,((m_1-m_2)/2.0)*4.359745e-6))
        plt.show()

    def plot_Heats(self,x2,datos):
        plt.plot(x2,self.Heat_1,x2,self.Heat_3)
        plt.xlabel('picoseconds')
        plt.ylabel('Hartree')

class VDOS_calculus(object):

    def read_xyz(self,filename_in,datos):
        """Lee las coordenadas de un archivo XYZ. Incluso las velocidades de acuerdo con el fichero geo_out.gen.
        Si el numero de coordenadas no coincide con lo establecido en el fichero da un aviso de error.
        Diego Martinez Gutierrez"""
        self.xyz = open(filename_in,"r")
        self.num_atoms = int(xyz.readline().split()[0])
        Vr = np.zeros((datos.num_pasos_tot,self.num_atoms)) # Vr[tiempo][atomo]
        Vr_X = np.zeros((datos.num_pasos_tot,self.num_atoms)) # Vr[tiempo][atomo]
        Vr_Y = np.zeros((datos.num_pasos_tot,self.num_atoms)) # Vr[tiempo][atomo]
        Vr_Z = np.zeros((datos.num_pasos_tot,self.num_atoms)) # Vr[tiempo][atomo]

        self.xyz_in = open(filename_in,"r")
        printProgress(0,datos.num_pasos_tot, prefix = 'Loading data:', suffix = ' ', barLength = 50)

        for tiempo in range(datos.num_pasos_tot):
            if (float(tiempo)/1000.0) == float(tiempo/1000): #progress bar no es necesario con pocos datos
                printProgress(tiempo,datos.num_pasos_tot, prefix = 'Loading data:', suffix = ' ', barLength = 50)
            num_atoms = int(self.xyz_in.readline().split()[0])
            title = self.xyz_in.readline()
            MD_step = title.split()[2]
            for atomo in range(num_atoms):
                atom,x,y,z,vx,vy,vz = self.xyz_in.readline().split()
                Vr[tiempo][atomo]=(float(vx)**2+float(vy)**2+float(vz)**2)**(0.5)
                Vr_X[tiempo][atomo]=float(vx)
                Vr_Y[tiempo][atomo]=float(vy)
                Vr_Z[tiempo][atomo]=float(vz)
        printProgress(0,self.num_atoms, prefix = 'Calculating VDOS:', suffix = ' ', barLength = 50)
        for atomo in datos.indice:
            if (float(atomo)/5.0) == float(atomo/5): #progress bar no es necesario con pocos datos
                printProgress(atomo,num_atoms, prefix = 'Calculating VDOS:', suffix = ' ', barLength = 50)
            if atomo == datos.indice[0]:
                VDOS = self.autocorrelation(Vr[datos.num_pasos_tot/2:,atomo])
                VDOS_X = self.autocorrelation_DOS(Vr_X[-datos.lenghtinterval::1,atomo])
                VDOS_Y = self.autocorrelation_DOS(Vr_Y[-datos.lenghtinterval::1,atomo])
                VDOS_Z = self.autocorrelation_DOS(Vr_Z[-datos.lenghtinterval::1,atomo])

            else:
                VDOS += self.autocorrelation(Vr[datos.num_pasos_tot/2:,atomo])
                VDOS_X += self.autocorrelation_DOS(Vr_X[-datos.lenghtinterval::1,atomo])
                VDOS_Y += self.autocorrelation_DOS(Vr_Y[-datos.lenghtinterval::1,atomo])
                VDOS_Z += self.autocorrelation_DOS(Vr_Z[-datos.lenghtinterval::1,atomo])

        xyz_in.close()

        return Vr,VDOS,VDOS_X,VDOS_Y,VDOS_Z
    def autocorrelation (self,x) :
        """
        Compute the autocorrelation of the signal, based on the properties of the
        power spectral density of the signal.
        """
        xp = x-np.mean(x)
        f = np.fft.fft(xp)
        p = np.array([v*np.conjugate(v) for v in f])
        pi = np.fft.ifft(p)
        return np.real(pi)[:x.size/2]/np.sum(xp**2)

    def autocorrelation_DOS(self,Vr):
        DOS = np.correlate(Vr,Vr,mode='full')
        return DOS

    def Plotter_VDOS(self,VDOS,VDOS_X,VDOS_Y,VDOS_Z,datos):
        plt.plot([x*33.35640/((datos.lenghtinterval*2.0)*datos.timestep) for x in range(len(VDOS_X))],((np.fft.fft(VDOS_X)*np.conjugate(np.fft.fft(VDOS_X))+np.fft.fft(VDOS_Y)*np.conjugate(np.fft.fft(VDOS_Y))+np.fft.fft(VDOS_Z)*np.conjugate(np.fft.fft(VDOS_Z)))**(0.5)))
        plt.xlabel('$\omega$ ($cm^{-1}$)')
        plt.ylabel('VDOS')
        plt.show()

        return 0

def main():
    datos = data()
#---------------------------part1:dT/dx---------------------------------
    parte1 = dTdx()
    parte1.read_xyz("geo_end.xyz","temp.file",datos)
    parte1.printGraphic(datos)
#--------------------------END:part1------------------------------------

#--------------------------Part2:dQ/dt----------------------------------
    parte2 = dQdt()
    parte2.read_md("md.out",parte2.Heat1,parte2.Heat2,parte2.Heat3,datos)
#--------------------------END:part2------------------------------------
    print('----------------------------------------------')
    print('G = {} W*m^-1*K^-1'.format(parte2.dQdt/(datos.Area*parte1.dTdx*1e-20)))
    print('G_k = {} W*m^-2*K^-1'.format(parte2.dQdt/(datos.Area*parte1.DeltaT)))
    print('G_k = {} GW*m^-2*K^-1'.format(parte2.dQdt*1e-9/(datos.Area*parte1.DeltaT)))
    print('Conductance(single chain) = {} W*K^-1'.format(parte2.dQdt/(parte1.DeltaT*datos.Cadenas)))
    print('Conductance(single chain) = {} nW*K^-1'.format(parte2.dQdt*1e9/(parte1.DeltaT*datos.Cadenas)))
    print('----------------------------------------------')
#--------------------------Part3:VDOS-----------------------------------
    parte3 = VDOS_calculus()

    Vr,VDOS,VDOS_X,VDOS_Y,VDOS_Z = parte3.read_xyz("geo_end.xyz",datos)
    parte3.Plotter_VDOS(VDOS,VDOS_X,VDOS_Y,VDOS_Z,datos)
#--------------------------END:Part3------------------------------------
    return 0

if __name__ == '__main__':
    main()
