#!/usr/bin/env python
# -*- coding: utf-8 -*-
# calculo de PbTe with SSCHA by Diego

# ---------------------------
# Importación de los módulos
# ---------------------------
from __future__ import print_function
from __future__ import division

import numpy
import matplotlib
from matplotlib import pylab, mlab, pyplot
np = numpy
plt = pyplot

from IPython.core.pylabtools import figsize, getfigs

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
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax



# -----------
# Constantes
# -----------

# ------------------------------
# Clases y Funciones utilizadas
# ------------------------------
class Mi_espresso(object):
    def __init__(self):
        #la parte del quantumespresso
        #definimos lo pseudopotenciales (con un diccionario)
        pseudos = {"Pb": "Pb.upf","Te": "Te.upf"}

        #definimos los parametros usados por quantumespresso
        input_params = {"ecutwfc" : 60, #cutoff de la funcion de onda en el plano (plane-wave wave-function cutoff)
        		"ecutrho": 240, #cutoff de la densidad de la funcion de onda
        		"conv_thr": 1e-6, #convergencia de la autoconsistencia DFT
        		"pseudo_dir": "pseudo_espresso", #directorio del pseudopotenciales
        		"tprnfor": True, #queremos que imprima las fuerzas
        		"tstress": True #queremos que imprima el tensor de estress
        		}

        k_spacing = 0.2 #distancia minima de sampleado en la celda Brillouin
        espresso_calc = Espresso(input_data = input_params, pseudopotentials = pseudos, kspacing = k_spacing)

    def relaja(self):
        #ahora realizamos el calculo mediante QuatumEspresso (ASE Espresso)
        #establecemos el tipo de calculo (una relajacion del cristal)
        input_params["calculation"] = "vc-relax"
        #de nuevo llamamos al Espresso
        espresso_calc = Espresso(input_data = input_params, pseudopotentials = pseudos, kspacing = k_spacing)
        #y lo unimos a nuestro sistema
        PbTe_primitive.set_calculator(espresso_calc)

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
        	nq1 = 2
        	nq2 = 2
        	nq3 = 2

        	! compute also the effective charges and the dielectric tensor
        	epsil = .true.
        &end
        """
        #ahora escribimos el fichero de entrada para el quantumespresso
        with open("harmonic.phi", "w") as f:
        	f.write(ph_input)
        #y enviamos el comando al sistema operativo para que haga el calculo, primero creamos el comando a ejecurar
        cmd = "mpirun -np 2 ph.x -nppolll 2 -i harmonic.phi > harmonic.pho" #por que limitar el calculo a dos nucleos solo?
        #cmd = "mpirun -np 8 ph.x -nppolll 2 -i harmonic.phi > harmonic.pho" #probemos con (cuatro) ocho nucleos, que es "nppoll"
        res = os.system(cmd) #esto es lo que pone en marcha el calculo, mirar a ver cuantos nucleos usa...

    def calcula(self):
    	#ahora calcularemos las fuerzas atomicas, energias BO, y el tensor de estres manualmente

    	#primero guardamos el 'ensemble' antes de hacer los calculos
    	#creamos el directorio de datos para guardar todos los ensembles de los pasos de la iteraccion del SSCHA
    	ensemble.save("data_ensemble_manual", population = 1) #creo que esto guarda en el directorio "data..." la poblacion 1 con las 10 configuraciones que hemos generado

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
    	    ntyp = 2
    	    ibrav = 0
    	    ecutwfc = 40
    	    ecutrho = 160
    	&end
    	&electrons
    	    conv_thr = 1d-6
    	    !diagonalization = "cg"
    	&end

    	ATOMIC_SPECIES
    	Pb 207.2 Pb.upf
    	Te 127.6 Te.upf

    	K_POINTS automatic
    	1 1 1 0 0 0
    	""".format(ensemble.structures[0].N_atoms)
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

def main()
    calcula_espresso = Mi_espresso()

    #importamos la estructura del PbTe
    PbTe_atoms = ase.io.read("PbTe.cif")

    #podemos ver la estructura
    #view(PbTe_atoms)

    #dado que la estructura no es de las dimensiones deseables para el SSCHA (dado que para la visualizacion se han declarado un exceso de atomos que a fin de calculos resultan antieconomicos)
    #redefinimos la celda unidad mediante el "CellConstructor"
    #inicializamos el CellConstructor con una nueva estructura
    struct = CC.Structure.Structure()
    struct.generate_from_ase_atoms(PbTe_atoms)

    #definimos la nueva celda unidad
    new_cell = struct.unit_cell.copy()
    #recontruimos los vectores v'1=(1/2)(v1+v2) v'2=0.5(v1-v2) v'3=0.5*v1+0.5*v3
    new_cell[0,:] = .5 * struct.unit_cell[0,:] + .5 * struct.unit_cell[1,:]
    new_cell[1,:] = .5 * struct.unit_cell[0,:] - .5 * struct.unit_cell[1,:]
    new_cell[2,:] = .5 * struct.unit_cell[0,:] + .5 * struct.unit_cell[2,:]

    #aplicamos la nueva celda unidad a la estructura y eliminamos duplicados
    struct.unit_cell = new_cell
    struct.fix_coords_in_unit_cell()
    PbTe_primitive = struct.get_ase_atoms()
    #podemos ver la estructura
    #view(PbTe_primitive)

    calcula_espresso.relaja()

    #ahora usamos el CellConstructor para importar los resultados de quantoespresso
    #para ello tenemos que saber cuantas matrices dinamicas ha obtenido espresso (en este caso han de ser 3)
    #cargamos las matrices dinamicas. he de mirar como lo hace...
    harm_dyn = CC.Phonons.Phonons("harmonic_dyn", nqirr = 3)

    #ahora imprimimos las frecuancias de los phonones
    w_s, pols = harm_dyn.DiagonalizeSupercell() #creo que esto calcula los valores propios (aka frecuencias)

    #aparentemente CellConstructor tiene incorporado insitu un conversor de unidades de Ry/bohr**2 a cm-1
    print ("\n".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in w_s]))

    #ahora llega la parte divertida SSCHA en accion!
    #leemos los datos de la matriz dinamica (tambien podemos hacer dyn=harm_dyn)
    dyn = CC.Phonons.Phonons("harmonic_dyn", nqirr = 3)

    #aplicamos las constricciones de la regla de la suma y las simetrias
    dyn.Symmetrize()

    #no queremos frecuencias imaginarias, hagamos que todo sea real. Ya bo es armonica pero en un estado inicial adecuado para SSCHA
    dyn.ForcePositiveDefinite()

    #veamos esas frecuencias, primero un separador
    print ("\n")
    w_s, pols = dyn.DiagonalizeSupercell() #de nuevo calcolamos los valores propios (frecuencias)
    print ("\n".join(["{:.4f} cm-1".format(w * CC.Units.RY_TO_CM) for w in w_s])) #ya sabemos como va esto...

    #la aproximacion estocastica (Monte Carlo)
    #nos preparamos para el SSCHA a T=100K tomando la matriz densidad de la matriz dinamica
    ensemble = sscha.Ensemble.Ensemble(dyn, T0 = 100, supercell = dyn.GetSupercell())
    #ensemble = sscha.Ensemble.Ensemble(dyn, T0 = 2000, supercell = dyn.GetSupercell())
    while True:
    	#creamos 10 estructuras para hacer Monte Carlo
    	ensemble.generate(N = 10)
    	#veamos como cambia al usar 20 estructuras <---------------------
    	#ensemble.generate(N = 20)
    	#visualizamos lo que hemos hecho
    	#view(ensemble.structures[0].get_ase_atoms())

    	#A fin de minimizar necesitamos conocer los potenciales, esto se hace usando algun codigo externo
    	#se puede hacer de forma manual o automatica

        calcula_espresso.calcula()

        #finalmente regresamos a SSCHA y cargamos los datos recordando que ahora tenemos que indicar el numero de configuraciones
    	ensemble.load("data_ensemble_manual", population = 1, N = 10)
    	#ahora con 20 poblaciones <-----------------------------------------------------------------
    	#ensemble.load("data_ensemble_manual", population = 1, N = 20)



    	#la minimizacion SSCHA

    	#reseteamos los calculos si se ejecuta esta celda en multiples ocasiones
    	ensemble.update_weights(dyn, 100) #restaura la matriz densidad a T=100k
    	minimizer = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

    	#se ignora la minimizacion de la estructura
    	minimizer.minim_struct = False

    	#se preparan los parametros para la minimizacion de la matriz covariante
    	minimizer.min_step_dyn = 0.05 #un valor de 1 es correcto
    	#...
    	#...

    	#se establece el parametro del chequeo de Kong-Liu, 0.5 es un buen comienzo, pero se pueden poner parametros mas bajos para acelerar la convergencia
    	minimizer.kong_liu_ratio = 0.3

    	#se hace correr la minimizacion
    	minimizer.init()
    	minimizer.run()

    	minimizer.plot_results()

    	minimizer.finalize()

    	#ahora de forma manual preguntamos si queremos continuar o paramos (si, esto es una chapuza)
    	if str(input("Do you wish to continue:[y/n]")) == 'n':
    		break #rompemos el bucle; veamos si funciona.

    #y en principio ahora solo queda mirar los resultados (cut-paste del manual):
    print("The total free energy per unit cell is:", minimizer.get_free_energy(), " Ry")
    print("The total stress tensor is [Ry/bohr^3]:")
    print(minimizer.get_stress_tensor()[0])
    print("And the stochastic error on the stress tensor is:")
    print(minimizer.get_stress_tensor()[1])
    print("The stocastic error of the free energy instead, was:", minimizer.get_free_energy(return_error = True)[1], " Ry")
    # Draw the 3D structure of the final average atomic positions
    #view(minimizer.dyn.structure.get_ase_atoms())

    # We can save the dynamical matrix
    minimizer.dyn.save_qe("dyn_pop1_")

    # Print the frequencies before and after the minimization
    w_old, p_old = ensemble.dyn_0.DiagonalizeSupercell() # This is the representation of the density matrix used to generate the ensemble
    w_new, p_new = minimizer.dyn.DiagonalizeSupercell()

    # We can now print them
    print(" Old frequencies |  New frequencies")
    print("\n".join(["{:16.4f} | {:16.4f}  cm-1".format(w_old[i] * CC.Units.RY_TO_CM, w_new[i] * CC.Units.RY_TO_CM) for i in range(len(w_old))]))

    #Ya que estamos podemos buscar inestabilidades estructurales por medio del Hessiano.
    #en este caso el hessiano de la energia libre se calcula con:
    #dyn_hessian = ensemble.get_free_energy_hessian()
    #por defecto se calcula con la aproximacion bubble, pero podemos incluir explicitamente los terminos  Phi4 con:
    #dyn_hessian = ensemble.get_free_energy_hessian(include_v4 = True)
    #tambien podemos poner verbose = True para obtener más información, incluyendo "d3_realspace_sym.npy"
    #dyn_hessian = ensemble.get_free_energy_hessian(verbose = True)
    dyn_hessian = ensemble.get_free_energy_hessian(include_v4 = True, verbose = True)
    dyn_hessian.save_qe("free_energy_hessian")

if __name__ == "__main__":
    main()
