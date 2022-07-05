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
import itertools

# -------
# Clases
# -------
class Signos(object):
    """docstring for Signos."""

    def __init__(self, arg):
        super(Signos, self).__init__()
        self.arg = arg
        self.tabla = []
        for v in itertools.product([+1,-1], repeat=6):
            #print(numpy.asarray(v))
            self.tabla.append(np.asarray(v))
        return 0

    def cut():
        """we will try to remove the antithetic conbinations """
        for i in range(len(tabla)):
            for j in range(i,len(tabla)):
                if ((tabla[i]*(-1)==tabla[j]).all()):
                    tabla = np.delete(tabla,j)  #does this work??
        return 0
    def unicos():
        new_array = [tuple(row) for row in tabla]
        uniques = np.unique(new_array)
        return uniques
# ----------
# Funciones
# ----------

def main(args):
    Signos()
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
