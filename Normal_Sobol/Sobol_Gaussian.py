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
# -------
# Clases
# -------
class Random_generador(object):
    def __init__(self,size):
        self.size=size
        ###########################################################
    def sobol_normal(self):
        sampler = qmc.Sobol(d=2, scramble=False)
        size_sobol = int(np.log(self.size)/np.log(2))+1
        print("size = ",self.size," ;size Sobol = ",size_sobol)
        sample = sampler.random_base2(m=size_sobol)

        data1=[]
        while (len(data1)<self.size):
         data = sampler.random()
         v1=2.0*data[0][0]-1.0
         v2=2.0*data[0][1]-1.0
         Riq=v1**2+v2**2
         if(0 <= Riq <= 1):
             #print("Riq=",Riq)
             data3=np.sqrt(-2*np.log(Riq)/Riq)
             #print("data=",v1*data3)
             data1.append(v1*data3)

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
         if(0 <= Riq <= 1):
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

def main(args):
    calculos = Random_generador(100)
    calculos.sobol_normal()
    calculos.random_normal()
    calculos.random_numpy_normal()
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
