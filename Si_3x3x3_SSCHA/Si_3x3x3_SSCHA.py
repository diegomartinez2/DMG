#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SSCHA_exercise.py
#
#  Copyright 2022 Diego <diego@u038025>
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

# Import Matplotlib to plot
import numpy as np
import matplotlib.pyplot as plt
# ---------------------------
# Importación de los módulos
# ---------------------------
from __future__ import print_function
from __future__ import division

import matplotlib
from matplotlib import pylab, mlab

#from IPython.core.pylabtools import figsize, getfigs

from pylab import *
from numpy import *


import sys,os
#importamos el metodo de calculo de las fuerzas Born-Oppenhaimer (Quantum-espresso)
import ase
from ase.calculators.espresso import Espresso
from ase.visualize import view
#importamos los modulos base del SSCHA (cellconstructor)
import cellconstructor as CC
import cellconstructor.Structure
import cellconstructor.Phonons

#importamos el motor SSCHA
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.Utilities

class Mi_espresso(object):
    def __init__(self):
        #la parte del quantumespresso
        #definimos lo pseudopotenciales (con un diccionario)
        self.pseudos = {"Si": "Si.upf"}

        #definimos los parametros usados por quantumespresso
        self.input_params = {"ecutwfc" : 60, #cutoff de la funcion de onda en el plano (plane-wave wave-function cutoff)
        		"ecutrho": 240, #cutoff de la densidad de la funcion de onda
        		"conv_thr": 1e-6, #convergencia de la autoconsistencia DFT
        		"pseudo_dir": "pseudo_espresso", #directorio del pseudopotenciales
        		"tprnfor": True, #queremos que imprima las fuerzas
        		"tstress": True #queremos que imprima el tensor de estress
        		}

        self.k_spacing = 0.2 #distancia minima de sampleado en la celda Brillouin
        self.espresso_calc = Espresso(input_data = self.input_params, pseudopotentials = self.pseudos, kspacing = self.k_spacing)
        def relaja(self,PbTe_primitive):
            #ahora realizamos el calculo mediante QuatumEspresso (ASE Espresso)
            #establecemos el tipo de calculo (una relajacion del cristal)
            self.input_params["calculation"] = "vc-relax"
            #de nuevo llamamos al Espresso
            self.espresso_calc = Espresso(input_data = self.input_params, pseudopotentials = self.pseudos, kspacing = self.k_spacing)
            #y lo unimos a nuestro sistema
            PbTe_primitive.set_calculator(self.espresso_calc)

            #y relajamos el cristal
            equilibrium_energy = PbTe_primitive.get_total_energy()
            #esto lleva un rato, nota porque en este ordenador no se usa mas de un nucleo?
            #ahora tenemos un  input file espresso.pwi
            #y un  output espresso.pwo
            #ahora la estructura ha cambiado, leemos del fichero de salida
            PbTe_final = ase.io.read("espresso.pwo")
            #y comprobamos que efectivamente se ha actualizado la estructura
            print("Volumen antes de la optimizacion:", PbTe_primitive.get_volume(), " A**3")
            print("Volumen despues de la optimizacion:", PbTe_final.get_volume(), " A**3")

            #harmonic calculation <-se necesita hacer un calculo armonico para el primer paso del SSCHA
            #creamos un input file para el quantumespresso para que haga el calculo de fonones
            ph_input = """
            &inputph
            	! the final filename
            	fildyn = "harmonic_dyn"

            	! the q mesh
            	ldisp = .true.
            	nq1 = 1
            	nq2 = 1
            	nq3 = 1

            	! compute also the effective charges and the dielectric tensor
            	epsil = .true.
            &end
            """
            #ahora escribimos el fichero de entrada para el quantumespresso
            with open("harmonic.phi", "w") as f:
            	f.write(ph_input)
            #y enviamos el comando al sistema operativo para que haga el calculo, primero creamos el comando a ejecurar
            #cmd = "mpirun -np 2 ph.x -nppolll 2 -i harmonic.phi > harmonic.pho" #por que limitar el calculo a dos nucleos solo?
            cmd = "mpirun -np 8 ph.x -nppolll 2 -i harmonic.phi > harmonic.pho" #probemos con (cuatro) ocho nucleos, que es "nppoll"
            res = os.system(cmd) #esto es lo que pone en marcha el calculo, mirar a ver cuantos nucleos usa...


def main():
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
