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


def File_Sort_Function(filename_sp,smear_cm, nsm):
    for ism in range(nsm):
        name = "{:6.2f}".format(smear_cm[ism]).strip()
        filename_data = filename_sp+'_'+name+'.dat'
        f = open(filename_data, 'rw')
        data = np.loadtxt(f)
        X = data[:,0]
        Y = data[:,1]
        data_out = [data[i] for i in np.lexsort((Y,X))] # <= Here is the sorting part
        np.savetxt(f,data_out)
    pass

def main(args):
    File_Sort_Function("test", 1.0)
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
