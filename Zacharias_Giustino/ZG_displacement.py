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

# -------
# Clases
# -------
class ZG_displacement(object):
    """Zacharias-Giustino displacement steps:
        1) Computing the vibrational eigenmodes and eigenfrequencies.
        2) Interpolation on q-points grid of eigenmodes and eigenfrequencies.
        3) Set the Berry connection.
        4) Calculate KG-displacement with known formula."""

    def __init__(self, dyn):
        super(ZG_displacement, self).__init__()
        self.new_dynamical_matrix = dyn.Copy()

    def ZG_eigen(self,new_dynamical_matrix):
        w, pols = self.new_dynamical_matrix.DiagonalizeSupercell()
        return w, pols

    def Interpolation(self, mesh_dim):
        self.dyn = self.new_dynamical_matrix.InterpolateMesh(mesh_dim)
        pass

    def Berry_connection(self,arg):
        pass

    def Signs(self, arg):
        pass

    def ZG_displacement(self, arg):
        """
        \Delta Tau_{p,k}=\sum_{q\in \mathcal{B},\nu} S_{q,\nu}[\frac{\hbar}{2N_p M_k w_{q\nu}}(2n_{q\nu,T}+1)]^{1/2}
        \times 2 Re[e^{iq\cdot R_p} pols_{k,\nu}(q)]
        """

        return Tau

# ----------
# Funciones
# ----------
def NombredeFuncion(arg):
    pass

def main(args):
    ZG_displace = ZG_displacement(dyn)
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
