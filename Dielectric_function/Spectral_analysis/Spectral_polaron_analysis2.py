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
import matplotlib.pyplot as plt
import numpy as np

class Polaron_analysis(object):
    """docstring fo Polaron_analysis."""

    def __init__(self, arg):
        super(Polaron_analysis, self).__init__()
        self.arg = arg

    def load_data(arg,namefile):
        with open(namefile) as file:
             lines = [line.rsplit() for line in file]
        Spec = np.zeros(255051)
        for i in range(255051):
             Spec[i]=lines[i][6]
        data = np.resize(Spec,(51,int(len(Spec)/51)))
        return data

    def plot_contour(arg,data):
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

def main(arg):
    namefile = "xbw"
    print ("Creando el objeto polaron")
    polaron = Polaron_analysis(arg)
    print ("Leyendo datos")
    data = polaron.load_data(namefile)
    print ("dibuja")
    polaron.plot_contour(data)
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
