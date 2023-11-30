#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Spectral_plasmons_analysis.py
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
import Thermal
import Spectral_function_calc

from __future__ import print_function
from __future__ import division
# Import the modules to read the dynamical matrix
import cellconstructor as CC
import cellconstructor.Phonons
# Import the SCHA modules
import sscha, sscha.Ensemble

import numpy as np
import cellconstructor.ForceTensor
import cellconstructor.ThermalConductivity
import time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def main(args):
    """
    """
    if (len( sys.argv ) > 1):
        N_RANDOM = arg[5]

        SAVE_PREFIX = 'dyn_hessian_'
        NQIRR = arg[2]
        Tg = arg[3]
        T =  arg[4]
        POPULATION = arg[1]
        DATA_DIR = "pop{}/data".format(POPULATION)
        DYN_PREFIX =  "pop{}/dyn_start_population3_".format(POPULATION)
        FINAL_DYN =   "pop{}/dyn_end_population3_".format(POPULATION)
        INCLUDE_V4 = False
        print ("Population:",POPULATION)
        print ("Number elements in the ensamble",N_RANDOM)
        print ("The number or irredubcible q points (nqirr):",NQIRR)
        print ("The temperature used to generate the configurations:",Tg)
        print ("The temperature for the calculation:",T)
        print ("Path to the directory ens_pop#lastpop where the last population is stored:",DATA_DIR)
        print ("The dynamical matrix that generated the last population:",DYN_PREFIX)
        print ("The SSCHA dynamical matrix obtained with the last minimization",FINAL_DYN)
        print ("Free energy Hessian dynamical matrices output:",SAVE_PREFIX)
        d3 = Hessian_calculus(DATA_DIR,N_RANDOM,DYN_PREFIX,FINAL_DYN,SAVE_PREFIX,
                        NQIRR,Tg,T,POPULATION,INCLUDE_V4)
        #-----
        Espectro =  Funcion_espectral(FINAL_DYN,NQIRR)
        Espectro.prepara_tensor()
        Espectro.calcula_espectro_basico_SrTiO3_multiprocessing(T,8)
        Espectro.calcula_full_correction_en_punto_G(T)
        Espectro.calcula_full_correction_en_punto_R(T)
        Espectro.calcula_oneshot_correction_en_punto_Gamma(T)
        Espectro.calcula_oneshot_correction_en_punto_R(T)
        Espectro.calcula_oneshot_correction_along_PATH(T)
        Espectro.calcula_espectro_correction_multiprocessing_SrTiO3(T,8)
        Espectro.dibuja_espectro_basico_SrTiO3(filename = "SrTiO3_static.dat", PATH = "GXMGRX")
        #-----
        thermal_calculo(d3,FINAL_DYN,NQIRR)
        harm_dos, anharm_dos = processing()
        plot(harm_dos, anharm_dos)
        np.savetxt("dos_harmonic.dat",harm_dos,header='Temperature dependent Harmonic DOS from auxiliary force constants:')
        np.savetxt("dos_anharmonic.dat",anharm_dos,header='Temperature dependent Anharmonic DOS from lineshapes: 2 lines of raw data 2 lines of gaussian smoothed data')

    else:
        print ("Arguments are Population, The number or irredubcible q points (nqirr), The temperature used to generate the configurations, The temperature for the calculation, and Number elements in the ensamble. In that order, separated by simple spaces")
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
