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
        super(NombredeClase, self).__init__()
        self.arg = arg

        self.amp1 = 50
        self.cen1 = 100
        self.wid1 = 5

        self.amp2 = 100
        self.cen2 = 150
        self.wid2 = 10

        self.amp3 = 50
        self.cen3 = 200
        self.wid3 = 5
        return 0

    def read_data(self):
        return 0

    def drawn_3D(self):
        return 0

    def drawn_2D(self):
        return 0

    def remove_noise(self):
        return 0

    def locate_3Lorenztians(self):
        return 0

# ----------
# Funciones
# ----------
def NombredeFuncion(arg):
    pass

def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
