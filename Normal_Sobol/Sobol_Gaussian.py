#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
#
#  Copyright 2022 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
import numpy as np
from scipy.stats import qmc
import matplotlib.pyplot as plt
import itertools, math

# -------
# Clases
# -------
class Random_generador(object):
    def __init__(self,size):
        self.size=size
        ###########################################################
    def sobol_normal(self,scramble=False):
        sampler = qmc.Sobol(d=2, scramble=scramble)
        size_sobol = int(np.log(self.size)/np.log(2))+1
        print("size = ",self.size," ;size Sobol = ",size_sobol)
        sample = sampler.random_base2(m=size_sobol)

        data1=[]
        while (len(data1)<self.size):
         data = sampler.random()
         v1=2.0*data[0][0]-1.0
         v2=2.0*data[0][1]-1.0
         Riq=v1**2+v2**2
#         if(0 <= Riq <= 1):
         if(0 < Riq <= 1):
             #print("Riq=",Riq)
             data3=np.sqrt(-2*np.log(Riq)/Riq)
             #print("data=",v1*data3)
             data1.append(v1*data3)
#         else:
#             print('Riq=',Riq)

        plt.hist(data1, bins=50)
        print("Sobol graphics")
        plt.savefig('Sobol_hist_{}.png'.format(self.size))
        plt.show()
        plt.scatter(data1,range(len(data1)))
        print("Sobol graphics")
        plt.savefig('Sobol_scatter_{}.png'.format(self.size))
        plt.show()
        return data1
        ###########################################################
    def random_normal(self):
        data1=[]
        while (len(data1)<self.size):
         v1=2.0*np.random.rand()-1.0
         v2=2.0*np.random.rand()-1.0
         Riq=v1**2+v2**2
         if(0 < Riq <= 1):
             #print("Riq=",Riq)
             data3=np.sqrt(-2*np.log(Riq)/Riq)
             #print("data=",v1*data3)
             data1.append(v1*data3)

        plt.hist(data1, bins=50)
        print("random normal graphics")
        plt.savefig('random_hist_{}.png'.format(self.size))
        plt.show()
        plt.scatter(data1,range(len(data1)))
        print("random normal graphics")
        plt.savefig('random_scatter_{}.png'.format(self.size))
        plt.show()
        return data1
        ###########################################################
    def random_numpy_normal(self):
        data1=[]
        while (len(data1)<self.size):
             data1.append(np.random.normal())
        plt.hist(data1, bins=50)
        print("numpy normal graphics")
        plt.savefig('numpy_normal_hist_{}.png'.format(self.size))
        plt.show()
        plt.scatter(data1,range(len(data1)))
        print("numpy normal graphics")
        plt.savefig('numpy_normal_scatter_{}.png'.format(self.size))
        plt.show()
        return data1
        ###########################################################
    def sobol_NEW_normal(self,scramble=False):
        sampler = qmc.Sobol(d=2, scramble=scramble)
        size_sobol = int(math.ceil(np.log(self.size/2)/np.log(2))) #int(np.log(self.size)/np.log(2))
        print("size = ",self.size," ;size Sobol = ",size_sobol)
        sample = sampler.random_base2(m=size_sobol)


        data1=[0,0]
        for i in (range(1,len(sample))):
            u1 = (sample[i][0]+0.01)%1 #****Diegom_test****
            u2 = (sample[i][1]+0.01)%1 # adding random number and modulus to avoid zeros
            r = -np.sqrt(-2*np.log(u1))
            theta = 2*np.pi*u2
            data1.append(r*np.cos(theta))
            data1.append(r*np.sin(theta))

        for i in (range(1,len(sample))):
            plt.scatter(sample[i][0],sample[i][1])
        plt.show()
        plt.hist(data1, bins=50)
        print("Sobol graphics:size=",self.size,'size_sobol=',size_sobol,'len data1=',len(data1))
        plt.savefig('Sobol_NEW_hist_{}.png'.format(self.size))
        plt.show()
        plt.scatter(data1,range(len(data1)))
        print("Sobol graphics")
        plt.savefig('Sobol_NEW_scatter_{}.png'.format(self.size))
        plt.show()
        print (np.resize(data1,(len(data1),3)))
        plt.hist(np.resize(data1,(len(data1),3)), bins=50)
        plt.show()
        return data1
# -------
# Functions
# -------
def sobol_norm_rand(size,n_modes,scramble=True): #***Diegom_test*** Not works, not a good distribution.
    sampler = qmc.Sobol(d=2,scramble=scramble)
    size_rand = size+n_modes
    size_sobol1 = int(np.log(size_rand)/np.log(2))+1
    sample = sampler.random_base2(m=size_sobol1)
    data1 = []
    while (len(data1<size_rand)):
        data = sampler.random()
        v1 = 2.0*data[0][0]-1.0
        v2 = 2.0*data[0][1]-1.0
        Riq = v1*v1+v2*v2
        if (0< Riq <= 1):
            data3 = np.sqrt(-2.0*np.log(Riq)/Riq)
            data1.append(v1*data3)
    return np.resize(data1,(size,n_modes))
################################################################################
def sobol_norm_rand2(size,n_modes,scramble=True):  # **** Diegom_test **** OK??
    sampler = qmc.Sobol(d=2,scramble=scramble)
    data2 = []
    size_sobol = int(np.log(size)/np.log(2))+1
    sample = sampler.random_base2(m=size_sobol)
    for i in range(n_modes):
        data1 = []
        while (len(data1)<size):
            data = sampler.random()  #+0.001*(np.random.randn()-1)
            v1 = 2.0*data[0][0]-1.0
            v2 = 2.0*data[0][1]-1.0
            Riq = v1*v1+v2*v2
            if (0< Riq <= 1):
                data3 = np.sqrt(-2.0*np.log(Riq)/Riq)
                data1.append(v1*data3)
        data2.append(data1)

    #return np.resize(data1,(size,n_modes))
    data4 = np.resize(data2,(size,n_modes))
    #plt.hist(data4[:,0], bins=50)
    for i in range(n_modes):
        plt.hist(data4[:,i], bins=50)
    #plt.hist(data1, bins=50)
    print("Sobol SSCHA graphics")
    plt.savefig('New_Sobol_hist_{}.png'.format(size))
    plt.show()
    plt.scatter(data4[:,0],range(len(data4[:,0])))
    #plt.scatter(data1,range(len(data1)))
    print("Sobol SSCHA graphics")
    plt.savefig('New_Sobol_scatter_{}.png'.format(size))
    plt.show()
    return np.resize(data1,(size,n_modes))
################################################################################
def sobol_norm_rand3(size,n_modes,scramble=False):  # **** Diegom_test **** do not work nicely.
        sampler = qmc.Sobol(d=2, scramble=scramble)
        size_sobol = int(np.log(size)/np.log(2))+1
        print("size = ",size," ;size Sobol = ",size_sobol)
        sample = sampler.random_base2(m=size_sobol)

        data1=[]
        i=0
        while (len(data1)<size):

         v1=2.0*sample[i][0]-1.0
         v2=2.0*sample[i][1]-1.0
         i += 1
         Riq=v1**2+v2**2
         if(0 < Riq <= 1):
             data3=np.sqrt(-2*np.log(Riq)/Riq)
             data1.append(v1*data3)
        data = data1*n_modes
        data1 = np.resize(data,(size,n_modes))
        for i in range(n_modes):
            plt.hist(data1[:,i], bins=50)
        print("Sobol3 graphics")
        plt.savefig('Sobol3_hist_{}.png'.format(size))
        plt.show()
        plt.scatter(data1[:,i],range(len(data1)))
        print("Sobol3 graphics")
        plt.savefig('Sobol3_scatter_{}.png'.format(size))
        plt.show()
        return data1
################################################################################
def sobol_norm_rand4(size,n_modes,scramble=True):  # **** Diegom_test **** OK??
    sampler = qmc.Sobol(d=2,scramble=scramble)
    data2 = []
    # size_sobol = int(np.log(size)/np.log(2))+1
    # sample = sampler.random_base2(m=size_sobol)
    for i in range(n_modes):
        data1 = []
        while (len(data1)<size):
            data = sampler.random()  #+0.001*(np.random.randn()-1)
            v1 = 2.0*data[0][0]-1.0
            v2 = 2.0*data[0][1]-1.0
            Riq = v1*v1+v2*v2
            if (0< Riq <= 1):
                data3 = np.sqrt(-2.0*np.log(Riq)/Riq)
                data1.append(v1*data3)
        data2.append(data1)

    #return np.resize(data1,(size,n_modes))
    data4 = np.resize(data2,(size,n_modes))
    #plt.hist(data4[:,0], bins=50)
    for i in range(n_modes):
        plt.hist(data4[:,i], bins=50)
    #plt.hist(data1, bins=50)
    print("Sobol SSCHA graphics")
    plt.savefig('New_Sobol4_hist_{}.png'.format(size))
    plt.show()
    plt.scatter(data4[:,0],range(len(data4[:,0])))
    #plt.scatter(data1,range(len(data1)))
    print("Sobol SSCHA graphics")
    plt.savefig('New_Sobol4_scatter_{}.png'.format(size))
    plt.show()
    return np.resize(data1,(size,n_modes))
################################################################################
def sobol_norm_rand5(size,n_modes,scramble=True):  # **** Diegom_test **** OK??
    # sampler = qmc.Sobol(d=2,scramble=scramble)
    # data2 = []
    # size_sobol = int(np.log(size)/np.log(2))+1
    # sample = sampler.random_base2(m=size_sobol)
    # data1 = []
    # for i in sample[:][0]:
    #     #data = sampler.random()  #+0.001*(np.random.randn()-1)
    #     v1 = 2.0*sample[i][0]-1.0
    #     v2 = 2.0*sample[i][1]-1.0
    #     Riq = v1*v1+v2*v2
    #     if (0< Riq <= 1):
    #         data3 = np.sqrt(-2.0*np.log(Riq)/Riq)
    #         data1.append(v1*data3)
    # data2.append(data1)
    #
    # #return np.resize(data1,(size,n_modes))
    # data4 = np.resize(data2,(size,n_modes))
    # #plt.hist(data4[:,0], bins=50)
    # for i in range(n_modes):
    #     plt.hist(data4[:,i], bins=50)
    # #plt.hist(data1, bins=50)
    # print("Sobol SSCHA graphics")
    # plt.savefig('New_Sobol5_hist_{}.png'.format(size))
    # plt.show()
    # plt.scatter(data4[:,0],range(len(data4[:,0])))
    # #plt.scatter(data1,range(len(data1)))
    # print("Sobol SSCHA graphics")
    # plt.savefig('New_Sobol5_scatter_{}.png'.format(size))
    # plt.show()
    # return np.resize(data1,(size,n_modes))
            sampler = qmc.Sobol(d=2, scramble=scramble)
            size_sobol = int(np.log(size)/np.log(2))+1
            print("size = ",size," ;size Sobol = ",size_sobol)
            sample = sampler.random_base2(m=size_sobol)

            data1=[]
            i=0
            for i in range(2**size_sobol):

             v1=2.0*sample[i][0]-1.0
             v2=2.0*sample[i][1]-1.0
             i += 1
             Riq=v1**2+v2**2
             if(0 < Riq <= 1):
                 data3=np.sqrt(-2*np.log(Riq)/Riq)
                 data1.append(v1*data3)
            plt.hist(data1, bins=50)
            plt.show()
            data = data1*n_modes
            data1 = np.resize(data,(size,n_modes))
            for i in range(n_modes):
                plt.hist(data1[:,i], bins=50)
            print("Sobol3 graphics")
            plt.savefig('Sobol3_hist_{}.png'.format(size))
            plt.show()
            plt.scatter(data1[:,i],range(len(data1)))
            print("Sobol3 graphics")
            plt.savefig('Sobol3_scatter_{}.png'.format(size))
            plt.show()
            return data1

def sobol_norm_rand5(size,n_modes,scramble=False):  # **** Diegom_test ****
    sampler = qmc.Sobol(d=2,scramble=scramble)
    size_sobol = int(np.log(size/2)/np.log(2))#+1
    size_new = 2**size_sobol
    sample = sampler.random_base2(m=size_sobol)
    random1 = [0]
    random2 = [0]
    for i in range(1,size_new):
        #print(i,sample[i][0],sample[i][1])
        r = -np.sqrt(-2.0 * np.log(sample[i][0]))         #<=== usa esto para crear la gaussiana
        theta = 2.0 * np.pi * sample[i][1]
        print ('i=',i,'r=',r,'theta=',theta,'random1=',r * np.cos(theta),'random2=',r * np.sin(theta))
        random1.append(r * np.cos(theta))
        random2.append(r * np.sin(theta))
    #print (random1)
    #print (random2)
    #random = (random1+random2)*n_modes
    random = np.resize((random1+random2)*n_modes,(size_new,n_modes))
    plt.hist(random1, bins=50)
    plt.show()
    plt.hist(random2, bins=50)
    plt.show()
    plt.hist(random, bins=50)
    plt.show()
    print ('size=',size,'size_new=',size_new,'return=',len(random1+random2))
    #plt.show()
    return random1+random2
################################################################################
def sobol_norm_rand6(size,n_modes,scramble=False):  # **** Diegom_test ****
    sampler = qmc.Sobol(d=2,scramble=scramble)
    size_sobol = int(math.ceil(np.log(size)/np.log(2)))#+1
    size_new = 2**size_sobol
    sample = sampler.random_base2(m=size_sobol)
    random = [0]
    for i in range(1,size_new,2):
        #print(i,sample[i][0],sample[i][1])
        r = -np.sqrt(-2.0 * np.log(sample[i][0]))         #<=== usa esto para crear la gaussiana
        theta = 2.0 * np.pi * sample[i][1]
        print ('i=',i,'r=',r,'theta=',theta,'random1=',r * np.cos(theta),'random2=',r * np.sin(theta))
        random.append(r * np.cos(theta))
        random.append(r * np.sin(theta))
    #print (random1)
    print (random,len(random))
    #random = (random1+random2)*n_modes
    random = np.resize(random,(size_new,n_modes))
    plt.hist(random, bins=50)
    plt.show()
    print ('size=',size,'size_new=',size_new,'return=',len(random))
    #plt.show()
    return random
################################################################################
def sobol_norm_rand7(size,n_modes,scramble=False):  # **** Diegom_test ****
        sampler = qmc.Sobol(d=2, scramble=scramble)
        size_sobol = int(math.ceil(np.log(size/2)/np.log(2))) #int(np.log(self.size)/np.log(2))
        print("size = ",size," ;size Sobol = ",size_sobol)
        sample = sampler.random_base2(m=size_sobol)


        data1=[0,0]
        for i in (range(1,len(sample))):
            u1 = sample[i][0]
            u2 = sample[i][1]
            r = -np.sqrt(-2*np.log(u1))
            theta = 2*np.pi*u2
            data1.append(r*np.cos(theta))
            data1.append(r*np.sin(theta))
        m = len(data1)-size
        data2 = data1[m:]
        print (len(data1),len(data2),size)
        data = np.resize(data2,(n_modes,size))

        for i in (range(1,len(sample))):
            plt.scatter(sample[i][0],sample[i][1])
        plt.show()
        for i in (range(n_modes)):
            plt.hist(data[i], bins=50)
        print("Sobol graphics:size=",size,'size_sobol=',size_sobol,'len data1=',len(data1))
        #plt.savefig('Sobol_NEW_hist_{}.png'.format(self.size))
        plt.show()
        for i in (range(n_modes)):
            plt.scatter(data[i],range(len(data[0])))
        print("Sobol graphics")
        #plt.savefig('Sobol_NEW_scatter_{}.png'.format(self.size))
        plt.show()
        return data
################################################################################
################################################################################
def main(args):
    # calculos = Random_generador(50)
    # calculos.sobol_normal(scramble=False)
    # calculos.sobol_NEW_normal(scramble=False) #!!!!!hay que terminarlo....
    # calculos.random_normal()
    # calculos.random_numpy_normal()
    sobol_norm_rand7(100,2,scramble=False)
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
