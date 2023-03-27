#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Average_variance_std.py
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
import numpy as np
# -------
# Clases
# -------
class NombredeClase(object):
    """docstring for NombredeClase."""

    def __init__(self, arg):
        super(NombredeClase, self).__init__()
        self.arg = arg

# ----------
# Funciones
# ----------
def ReadList(arg):
#    file1 = open("MyFile.txt", "r")
#    list=file1.readlines()
    separated_lines = [line.strip() for line in open(arg)]
#    flattened_list = [item for i in separated_lines for item in i.split(",")]
#    file1.close()
    list=[float(line) for line in separated_lines]
    return list

def main(args):
    list=ReadList(args[1])
    print("average= ",np.average(list))
    print("variance= ",np.var(list))
    print("standar deviation= ",p.std(list))
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
