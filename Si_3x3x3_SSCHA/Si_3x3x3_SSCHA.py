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
from __future__ import print_function
from __future__ import division
# ---------------------------
# Importación de los módulos
# ---------------------------
# Import Matplotlib to plot
import numpy as np
import matplotlib.pyplot as plt

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
        #self.pseudos = {"Si": "Si.upf"}
        self.pseudos = {"Si": "Si.pbe-n-rrkus_psl.1.0.0.UPF"}

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
    def relaja(self,Si_primitive):
        #ahora realizamos el calculo mediante QuatumEspresso (ASE Espresso)
        #establecemos el tipo de calculo (una relajacion del cristal)
        self.input_params["calculation"] = "vc-relax"
        #de nuevo llamamos al Espresso
        self.espresso_calc = Espresso(input_data = self.input_params, pseudopotentials = self.pseudos, kspacing = self.k_spacing)
        #y lo unimos a nuestro sistema
        Si_primitive.set_calculator(self.espresso_calc)

        #y relajamos el cristal
        equilibrium_energy = Si_primitive.get_total_energy()
        #esto lleva un rato, nota porque en este ordenador no se usa mas de un nucleo?
        #ahora tenemos un  input file espresso.pwi
        #y un  output espresso.pwo
        #ahora la estructura ha cambiado, leemos del fichero de salida
        Si_final = ase.io.read("espresso.pwo")
        #y comprobamos que efectivamente se ha actualizado la estructura
        print("Volumen antes de la optimizacion:", Si_primitive.get_volume(), " A**3")
        print("Volumen despues de la optimizacion:", Si_final.get_volume(), " A**3")

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

    def calculo(self,ensemble,poblacion):
        #ahora calcularemos las fuerzas atomicas, energias BO, y el tensor de estres manualmente
    	#primero guardamos el 'ensemble' antes de hacer los calculos
    	#creamos el directorio de datos para guardar todos los ensembles de los pasos de la iteraccion del SSCHA
        ensemble.save("data_ensemble_manual", population = poblacion) #creo que esto guarda en el directorio "data..." la poblacion 1 con las 10 configuraciones que hemos generado
    	#ahora usamos python-magik para crear el fichero de entrada del motor de calculo externo
    	#Header (no tengo mas que una idea basica de lo que pone aqui salvo "scf" es el metodo, nat es el numero de atomos, numero de tipos de atomos 2 (ntyp) ... )
        typical_espresso_header = """
        &control
    	    calculation = "scf"
    	    tstress = .true.
    	    tprnfor = .true.
    	    disk_io = "none"
    	    pseudo_dir = "pseudo_espresso"
    	&end
    	&system
    	    nat = {}
    	    ntyp = 1
    	    ibrav = 0
    	    ecutwfc = 40
    	    ecutrho = 160
    	&end
    	&electrons
    	    conv_thr = 1d-6
    	    !diagonalization = "cg"
    	&end

    	ATOMIC_SPECIES
    	Si 28,085 Si.pbe-n-rrkus_psl.1.0.0.UPF
    	K_POINTS automatic
    	1 1 1 0 0 0
    	""".format(ensemble.structures[0].N_atoms)
        #Si 28,085 Si.upf
        #'Si'    47.867     'Si.pbesol-n-rrkjus_psl.0.1.UPF'

    	#el numero de atomos lol extraemos de los datos internos del ensemble
    	# We extract the number of atoms form the ensemble and the celldm(1) from the dynamical matrix (it is stored in Angstrom, but espresso wants it in Bohr)
    	# You can also read it on the fourth value of the first data line on the first dynamical matrix file (dyn_start_popilation1_1); In the latter case, it will be already in Bohr.

    	#llemos los datos
        all_scf_files = [os.path.join("data_ensemble_manual", f) for f in os.listdir("data_ensemble_manual") if f.startswith("scf_")]

    	#generamos un input file en un nuevo directorio si no existe ya
        if not os.path.exists("run_calculation"):
            os.mkdir("run_calculation")
        for file in all_scf_files:
            number = int(file.split("_")[-1].split(".")[0]) #muy inteligente, extraemos el X de "data_ensemble_manual/scf_population1_X.dat" primero separamos en dos grupos separados por "." y escogemos la parte de delante [0] y luego separamos los grupos divididos en "_" y tomamops el ultimo valor [-1]
            filename = os.path.join("run_calculation", "espresso_run_{}.pwi".format(number)) #usamos el numero recien obtenido para crear el nobre del fichero de entrada de los calculos
            with open(filename, "w") as f: #aqui escribimos los ficheros de entrada
                f.write(typical_espresso_header) #nuestro headeer va aqui
                ff = open(file, "r")	#abrimos el fichero de datos
                structure_lines = ff.readlines()	#pasamos las lineas del fichero a una variable structure_lines
                ff.close() #cerramos fichero
                f.writelines(structure_lines) #añadimos los datos al fichero de entrada

    	#con esto tenemos los ficheros de entrada necesarios para hacer un calculo en remoto.
    	#---------hasta aquí nsi se quere hacer en un cluster lejano------
    	#pero tambien lo podemos hacer local:

        directory = "run_calculation"
        for file in os.listdir(directory):
            if not file.endswith(".pwi"): #si no es parte del imput del calculo
                continue	#pasamos de ello
            outputname = file.replace(".pwi",".pwo") #la salida tiene el mismo nobre pero con la extension cambiada
            total_inputname = os.path.join(directory, file) #le ponemos la ruta del directorio
            total_outputname = os.path.join(directory, outputname) #le ponemos la ruta del directorio
    		#ahora hacemos correr el calculo con (cuatro) 8 procesadores
    	#	cmd = "mpirun -np 4 pw.x -i {} > {}".format(total_inputname, total_outputname)
            cmd = "mpirun -np 8 pw.x -i {} > {}".format(total_inputname, total_outputname)
            print("Ahora corriendo: ", cmd) #comienza el calculo
            os.system(cmd) #y finalmente se ejecuta
    		#probando sin mpi
    	#	cmd2 = "pw.x -i {} > {}".format(total_inputname, total_outputname)
    	#	print("Ahora corriendo: ", cmd2) #comienza el calculo
    	#	os.system(cmd2)
    	#--------fin del calculo local------
    	#ahora por local o por correo tenemos en run_calculation los ficheros de salida del metodo de calculo de fuerzas, presion, etc.
        directory = "run_calculation"
        output_filenames = [f for f in os.listdir(directory) if f.endswith(".pwo")]
        output_files = [os.path.join(directory, f) for f in output_filenames] #si, todo esto solo para saber que ficheros son y donde estan

    	#preparamos el array de las energias
        energies = np.zeros(len(output_files)) #tantas energias como ficheros
        for file in output_files:
            id_number = int(file.split("_")[-1].split(".")[0]) #cargamos el numero de la configuracion
            ff = open(file, "r") #vamos a leer el fichero
            lines = [l.strip() for l in ff.readlines()] #leemos las lineas y de paso quitamos los espacios en blanco que sobran
            ff.close()

    		#en espresso la primera linea del archivo que tien un "!" es la energia
            energy_line = next(l for l in lines if len(l) > 0 if l.split()[0] == "!")
            energies[id_number -1] = float(energy_line.split()[4]) #el 5 elemento de la linea de energia es la energia y la metemos en el array (recuerda que en python el array empieza por 0)
    		#ahora el numero de atomos
    	#	nat_line = next(l for l in lines if len(l) > 0 if l.split()[0] == "number" and l.split()[2] == "atoms/cell" )
            nat_line = next( l for l in lines if len(l) > 0 if l.split()[0] == "number" and l.split()[2] == "atoms/cell" )
            nat = int(nat_line.split()[4]) #es tan complejo extraer datos de un fichero de texto con formato variable...
    		#ahora las fuerzas
            forces = np.zeros((nat, 3))
            forces_lines = [l for l in lines if len(l) >0 if l.split()[0] == "atom"]
            for i in range(nat):
                forces[i,:] = [float(x) for x in forces_lines[i].split()[-3:]] #las fuerzas son los tres ultimos numeros
    		#ahora el stress
            stress = np.zeros((3,3))
    		#como siempre el formato del fichero es diferente aqui tambien primero tiene una cabecera y luego los datos
            index_before_stress = next(i for i, l in enumerate(lines) if len(l) > 0 if l.split()[0] == "total" and l.split()[1] == "stress") #recuerdese que enumerate() da duplas de datos (i,j) y solo queremos en que linea esta
            for i in range(3):
                index = i + index_before_stress +1
                stress[i, :] = [float(x) for x in lines[index].split()[:3]]

    		#guardamos los datos porque extraerlos ha sido muy molesto y no quiero hacerlo otra vez
            force_file = os.path.join("data_ensemble_manual", "forces_population1_{}.dat".format(id_number))
            stress_file = os.path.join("data_ensemble_manual", "pressures_population1_{}.dat".format(id_number))
            np.savetxt(force_file, forces)
            np.savetxt(stress_file, stress)

    	#guardamos los datos de la energia ... finalmente
        energy_file = os.path.join("data_ensemble_manual", "energies_supercell_population1.dat")
        np.savetxt(energy_file, energies)


def main(args):
    calcula_espresso = Mi_espresso()
    Si_atoms = ase.io.read("Si.cif")
    struct = CC.Structure.Structure()
    struct.generate_from_ase_atoms(Si_atoms)
    struct.fix_coords_in_unit_cell()
    Si_primitive = struct.get_ase_atoms()
    calcula_espresso.relaja(Si_primitive)


    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
