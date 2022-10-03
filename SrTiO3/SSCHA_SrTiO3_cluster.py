#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SSCHA_SrTiO3_cluster.py
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
# Import the sscha code
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax, sscha.Utilities

# Import the cellconstructor library to manage phonons
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.Structure, cellconstructor.calculators

# Import the force field of Gold
import ase, ase.calculators
from ase.calculators.emt import EMT

# Import numerical and general pourpouse libraries
import numpy as np, matplotlib.pyplot as plt
import sys, os


import cellconstructor.ForceTensor
import ase.dft.kpoints

import matplotlib.cm as cm
import matplotlib.colors as colors

import scipy, scipy.optimize

import sscha.Cluster

class Send_to_cluster(object):
    def __init__(self,hostname = 'diegom@ekhi.cfm.ehu.es', pwd = None,
           label = "", account_name = '', n_nodes = 1,
           time = '12:00:00', n_pool = 1, workdir = "/scratch/diegom/SrTiO3"):
        self.cluster = sscha.Cluster.Cluster(hostname = hostname, pwd = pwd)  # Put the password in pwd if needed

        # Configure the submission strategy
        if (account_name != ''):
            self.cluster.account_name = account_name  # Name of the account on which to subtract nodes
        if (account_name == ''):
            self.cluster.use_account = False
        self.cluster.n_nodes = n_nodes            # Number of nodes requested for each job
        self.cluster.time = time                  # Total time requested for each job
        self.cluster.n_pool = n_pool              # Number of pools for the Quantum ESPRESSO calculation
        self.cluster.label = label
        # Here some custom parameters for the clusters
        # These are specific for daint, but you can easily figure out those for your machine
        #self.cluster.custom_params["--constraint"] = "gpu"      # Run on the GPU partition
        #self.cluster.custom_params["--ntasks-per-node"] = '2'
        #self.cluster.custom_params["--cpus-per-task"] = '6'

        # Since daint specify the partition with a custom option,
        # Lets remove the specific partition option of SLURM
        # Neither we want to specify the total number of cpus (automatically determined by the node)
        self.cluster.use_partition = False
        self.cluster.use_cpu = False


        # Now, we need to tell daint which modules to load to run quantum espresso
        # Also this is cluster specific, but very simple to figure it out for you
        self.cluster.load_modules = """
        # Load the quantum espresso modules
        ## module load daint-gpu
        module load QuantumESPRESSO

        # Configure the environmental variables of the job
        export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
        ## export NO_STOP_MESSAGE=1
        ## export CRAY_CUDA_MPS=1

        ## ulimit -s unlimited
        """

        # Now, what is the command to run quantum espresso on the cluster?
        self.cluster.binary = "pw.x -npool NPOOL -i PREFIX.pwi > PREFIX.pwo"  #<--- No need for mpirun??
        #self.cluster.binary = "mpirun -np %d pw.x -npool NPOOL -i PREFIX.pwi > PREFIX.pwo" % mpi_kernels
        #self.cluster.binary = "mpirun -np {np} pw.x -npool NPOOL -i PREFIX.pwi > PREFIX.pwo".format(np = 8)
        # NOTE that NPOOL will be replaced automatically with the cluster.n_pool variable


        # Let us setup the working directory (directory in which the jobs runs)
        #self.cluster.workdir =  "$SCRATCH/diegom/Gold_NVT_300k"       # <--- IOf this fails just use "/scratch/[...]" instead
        self.cluster.workdir = workdir
        self.cluster.setup_workdir()  # Login to the cluster and creates the working directory if it does not exist


        # Last but not least:
        #  How many jobs do you want to submit simultaneously?
        self.cluster.batch_size = 10

        #  Now many DFT calculation do you want to run inside each job?
        self.cluster.job_number = 1 #5 before

class Gold_free_energy(object):
    def __init__(self):
        """
        Here we load the primitive cell of Gold from a cif file.
        And we use CellConstructor to compute phonons from finite differences.
        The phonons are computed on a q-mesh 4x4x4
        """

        self.gold_structure = CC.Structure.Structure()
        self.gold_structure.read_generic_file("Au.cif")

        # Get the force field for gold
        self.calculator = EMT()
        #self.Calculadora = Dft_calculator()
        #self.calculator = self.Calculadora.calculator

    def relax(self):
        # Relax the gold structure (useless since for symmetries it is already relaxed)
        relax = CC.calculators.Relax(self.gold_structure, self.calculator)
        self.gold_structure_relaxed = relax.static_relax()
        #return 0

    def harmonic(self):
        # Compute the harmonic phonons
        # NOTE: if the code is run with mpirun, the calculation goes in parallel
        gold_harmonic_dyn = CC.Phonons.compute_phonons_finite_displacements(
              self.gold_structure_relaxed, self.calculator, supercell = (4,4,4))

        # Impose the symmetries and
        # save the dynamical matrix in the quantum espresso format
        gold_harmonic_dyn.Symmetrize()
        gold_harmonic_dyn.save_qe("harmonic_dyn")


        # If the dynamical matrix has imaginary frequencies, remove them
        gold_harmonic_dyn.ForcePositiveDefinite()

        """
        gold_harmonic_dyn is ready to start the SSCHA calculation.

        Now let us initialize the ensemble, and the calculation at 300 K.
        We will run a NVT calculation, using 100 configurations at each step
        """

        TEMPERATURE = 300
        N_CONFIGS = 50
        MAX_ITERATIONS = 20

        # Initialize the random ionic ensemble
        ensemble = sscha.Ensemble.Ensemble(gold_harmonic_dyn, TEMPERATURE)

        # Initialize the free energy minimizer
        minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
        minim.set_minimization_step(0.01)

        # Initialize the NVT simulation
        mi_cluster = Send_to_cluster()
        relax = sscha.Relax.SSCHA(minim, self.calculator, N_configs = N_CONFIGS,
                         max_pop = MAX_ITERATIONS, cluster = mi_cluster.cluster,
                         save_ensemble = True)

        # Define the I/O operations
        # To save info about the free energy minimization after each step
        ioinfo = sscha.Utilities.IOInfo()
        ioinfo.SetupSaving("minim_info")
        relax.setup_custom_functions(custom_function_post = ioinfo.CFP_SaveAll)


        # Run the NVT simulation (save the stress to compute the pressure)
        relax.relax(get_stress = True,sobol = True)

        # If instead you want to run a NPT simulation, use
        # The target pressure is given in GPa.
        #relax.vc_relax(target_press = 0)

        # You can also run a mixed simulation (NVT) but with variable lattice parameters
        #relax.vc_relax(fix_volume = True)

        # Now we can save the final dynamical matrix
        # And print in stdout the info about the minimization
        relax.minim.finalize()
        relax.minim.dyn.save_qe("sscha_T{}_dyn".format(TEMPERATURE))
        return 0
class  Gold_free_energy_ab_initio(object):
    def __init__(self):
        # Initialize the DFT (Quantum Espresso) calculator for gold
        # The input data is a dictionary that encodes the pw.x input file namelist
        input_data = {
            'control' : {
                # Avoid writing wavefunctions on the disk
                'disk_io' : 'None',
                # Where to find the pseudopotential
                #'pseudo_dir' : '.'
            },
            'system' : {
                # Specify the basis set cutoffs
                'ecutwfc' : 45,   # Cutoff for wavefunction
                'ecutrho' : 45*4, # Cutoff for the density
                # Information about smearing (it is a metal)
                'occupations' : 'smearing',
                'smearing' : 'mv',
                'degauss' : 0.03
            },
            'electrons' : {
                'conv_thr' : 1e-8
            }
        }

        # the pseudopotential for each chemical element
        # In this case just Gold
        pseudopotentials = {'Au' : 'Au_ONCV_PBE-1.0.oncvpsp.upf'}

        # the kpoints mesh and the offset
        kpts = (1,1,1)
        koffset = (1,1,1)

        # Specify the command to call quantum espresso
        command = 'mpirun -np 8 pw.x -i PREFIX.pwi > PREFIX.pwo'


        # Prepare the quantum espresso calculator
        self.calculator = CC.calculators.Espresso(input_data,
                                             pseudopotentials,
                                             command = command,
                                             kpts = kpts,
                                             koffset = koffset)

    def relax(self):

        TEMPERATURE = 300
        N_CONFIGS = 32
        MAX_ITERATIONS = 20
        START_DYN = 'start_dyn'
        NQIRR = 13

        # Let us load the starting dynamical matrix
        gold_dyn = CC.Phonons.Phonons(START_DYN, NQIRR)

        # Initialize the random ionic ensemble
        ensemble = sscha.Ensemble.Ensemble(gold_dyn, TEMPERATURE)

        # Initialize the free energy minimizer
        minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
        minim.set_minimization_step(0.01)

        # Initialize the NVT simulation
        mi_cluster = Send_to_cluster()
        relax = sscha.Relax.SSCHA(minim, self.calculator, N_configs = N_CONFIGS,
                         max_pop = MAX_ITERATIONS, cluster = mi_cluster.cluster,
                         save_ensemble = True)

        # Define the I/O operations
        # To save info about the free energy minimization after each step
        ioinfo = sscha.Utilities.IOInfo()
        ioinfo.SetupSaving("minim_info")
        relax.setup_custom_functions(custom_function_post = ioinfo.CFP_SaveAll)


        # Run the NVT simulation (save the stress to compute the pressure)
        relax.relax(get_stress = True, sobol = True)

        # If instead you want to run a NPT simulation, use
        # The target pressure is given in GPa.
        #relax.vc_relax(target_press = 0)

        # You can also run a mixed simulation (NVT) but with variable lattice parameters
        #relax.vc_relax(fix_volume = True)

        # Now we can save the final dynamical matrix
        # And print in stdout the info about the minimization
        relax.minim.finalize()
        relax.minim.dyn.save_qe("sscha_T{}_dyn".format(TEMPERATURE))

class  Thermal_expansion(object):
    def __init__(self):
        """
        You need first to run the
        get_gold_free_energy.py

        Here we use NPT simulation to compute the gold thermal expansion.
        """

        # Define the temperature range (in K)
        self.T_START = 300
        self.T_END = 1000
        self.DT = 50

        self.N_CONFIGS = 32
        self.MAX_ITERATIONS = 10

        # Import the gold force field
        self.calculator = EMT()
        #self.Calculadora = Dft_calculator()
        #self.calculator = self.Calculadora.calculator

        # Import the starting dynamical matrix (final result of get_gold_free_energy.py)
        self.dyn = CC.Phonons.Phonons("sscha_T300_dyn", nqirr = 13)

        # Create the directory on which to store the output
        self.DIRECTORY = "thermal_expansion"
        if not os.path.exists(self.DIRECTORY):
            os.makedirs("thermal_expansion")

    def temperature_cycle(self):
        # We cycle over several temperatures
        t = self.T_START


        volumes = []
        temperatures = []
        while t <= self.T_END:
            # Change the temperature
            ensemble = sscha.Ensemble.Ensemble(self.dyn, t)
            minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
            minim.set_minimization_step(0.1)
            mi_cluster = Send_to_cluster()
            relax = sscha.Relax.SSCHA(minim, self.calculator,
                                        N_configs = self.N_CONFIGS,
                                        max_pop = self.MAX_ITERATIONS,
                                        cluster = mi_cluster.cluster,
                                        save_ensemble = True)

            # Setup the I/O
            ioinfo = sscha.Utilities.IOInfo()
            ioinfo.SetupSaving( os.path.join(self.DIRECTORY, "minim_t{}".format(t)))
            relax.setup_custom_functions( custom_function_post = ioinfo.CFP_SaveAll)


            # Run the NPT simulation
            relax.vc_relax(target_press = 0, sobol = True)

            # Save the volume and temperature
            volumes.append(relax.minim.dyn.structure.get_volume())
            temperatures.append(t)

            # Start the next simulation from the results
            relax.minim.dyn.save_qe( os.path.join(self.DIRECTORY, "sscha_T{}_dyn".format(t)))
            dyn = relax.minim.dyn
            relax.minim.finalize()

            # Update the temperature
            t += self.DT

        # Save the thermal expansion
        np.savetxt(os.path.join(self.DIRECTORY, "thermal_expansion.dat"),
                   np.transpose([temperatures, volumes]),
                   header = "Temperature [K]; Volume [A^3]")
        #return 0

class Dft_calculator(object):
    def __init__(self):

        # Initialize the DFT (Quantum Espresso) calculator for gold
        # The input data is a dictionary that encodes the pw.x input file namelist
        input_data = {
            'control' : {
                # Avoid writing wavefunctions on the disk
                'disk_io' : 'None',
                # Where to find the pseudopotential
                #'pseudo_dir' : '.'
            },
            'system' : {
                # Specify the basis set cutoffs
                'ecutwfc' : 45,   # Cutoff for wavefunction
                'ecutrho' : 45*4, # Cutoff for the density
                # Information about smearing (it is a metal)
                'occupations' : 'smearing',
                'smearing' : 'mv',
                'degauss' : 0.03
            },
            'electrons' : {
                'conv_thr' : 1e-8
            }
        }

        # the pseudopotential for each chemical element
        # In this case just Gold
        pseudopotentials = {'Au' : 'Au_ONCV_PBE-1.0.oncvpsp.upf'}

        # the kpoints mesh and the offset
        kpts = (1,1,1)
        koffset = (1,1,1)

        # Specify the command to call quantum espresso
        self.command = 'pw.x -i PREFIX.pwi > PREFIX.pwo'

        # Prepare the quantum espresso calculator
        self.calculator = CC.calculators.Espresso(input_data,
                                             pseudopotentials,
                                             kpts = kpts,
                                             koffset = koffset)


def plot_dispersion():
    NQIRR = 13
    #CMAP = "Spectral_r"
    PATH = "GXWXKGL"
    N_POINTS = 1000

    SPECIAL_POINTS = {"G": [0,0,0],
                      "X": [0, .5, .5],
                      "L": [.5, .5, .5],
                      "W": [.25, .75, .5],
                      "K": [3/8., 3/4., 3/8.]}

    # Load the harmonic and sscha phonons
    harmonic_dyn = CC.Phonons.Phonons('harmonic_dyn', NQIRR)
    sscha_dyn = CC.Phonons.Phonons('sscha_T300_dyn', NQIRR)

    # Get the band path
    qpath, data = CC.Methods.get_bandpath(harmonic_dyn.structure.unit_cell,
                                          PATH,
                                          SPECIAL_POINTS,
                                          N_POINTS)
    xaxis, xticks, xlabels = data # Info to plot correclty the x axis

    # Get the phonon dispersion along the path
    harmonic_dispersion = CC.ForceTensor.get_phonons_in_qpath(harmonic_dyn, qpath)
    sscha_dispersion = CC.ForceTensor.get_phonons_in_qpath(sscha_dyn, qpath)

    nmodes = harmonic_dyn.structure.N_atoms * 3

    # Plot the two dispersions
    plt.figure(dpi = 150)
    ax = plt.gca()

    for i in range(nmodes):
        lbl=None
        lblsscha = None
        if i == 0:
            lbl = 'Harmonic'
            lblsscha = 'SSCHA'

        ax.plot(xaxis, harmonic_dispersion[:,i], color = 'k', ls = 'dashed', label = lbl)
        ax.plot(xaxis, sscha_dispersion[:,i], color = 'r', label = lblsscha)

    # Plot vertical lines for each high symmetry points
    for x in xticks:
        ax.axvline(x, 0, 1, color = "k", lw = 0.4)
    ax.axhline(0, 0, 1, color = 'k', ls = ':', lw = 0.4)

    ax.legend()

    # Set the x labels to the high symmetry points
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

    ax.set_xlabel("Q path")
    ax.set_ylabel("Phonons [cm-1]")

    plt.tight_layout()
    plt.savefig("dispersion.png")
    plt.show()
    return 0

def plot_thermal_expansion():
    """
    This simple scripts plot the thermal expansion.
    It either directly plots the file produced at the end of thermal_expansion.py
    or it processes the dynamical matrices generated throughout the minimization.

    It also fits the V(T) cuve and estimates the volumetric thermal expansion coefficient

    alpha = dV/dT / V

    At 300 K
    """


    # Load all the dynamical matrices and compute volume
    DIRECTORY = "thermal_expansion"
    FILE = os.path.join(DIRECTORY, "thermal_expansion.dat")


    # Check if the file with the thermal expansion data exists
    if not os.path.exists( FILE):
        # If the simulation is not ended, load the volume from the dynamical matrices

        # Get the dynamical matrices
        all_dyn_files = [x for x in os.listdir(DIRECTORY) if "sscha" in x and x.endswith("dyn1")]
        temperatures = [float(x.split("_")[-2].replace("T", "")) for x in all_dyn_files]

        # Now sort in order of temperature
        sortmask = np.argsort(temperatures)
        all_dyn_files = [all_dyn_files[x] for x in sortmask]
        temperatures = np.sort(temperatures)

        volumes = np.zeros_like(temperatures)

        for i, fname in enumerate(all_dyn_files):
            # Load the dynamical matrix
            # The full_name means that we specify the name including the tail 1
            dyn = CC.Phonons.Phonons(os.path.join(DIRECTORY, fname), full_name = True)
            volumes[i] = dyn.get_volumes()

    else:
        # Load the data from the final data file
        temperatures, volumes = np.loadtxt(FILE, unpack = True)


    # Prepare the figure and plot the V(T) from the sscha data
    plt.figure(dpi = 150)
    plt.scatter(temperatures, volumes, label = "SSCHA data", color = 'r')

    # Fit the data with a quadratic curve
    def parabola(x, a, b, c):
        return a + b*x + c*x**2
    def diff_parab(x, a, b, c):
        return b + 2*c*x

    popt, pcov = scipy.optimize.curve_fit(parabola, temperatures, volumes,
                                          p0 = [0,0,0])

    # Evaluate the volume thermal expansion
    vol_thermal_expansion = diff_parab(300, *popt) / parabola(300, *popt)
    print("Vol thermal expansion: {} x 10^6  K^-1".format(vol_thermal_expansion * 1e6))
    plt.text(0.6, 0.2, r"$\alpha_v = "+"{:.1f}".format(vol_thermal_expansion*1e6)+r"\times 10^6 $ K$^{-1}$",
             transform = plt.gca().transAxes)


    # Plot the fit
    _t_ = np.linspace(np.min(temperatures), np.max(temperatures), 1000)
    plt.plot(_t_, parabola(_t_, *popt), label = "Fit", color = 'k', zorder = 0)

    # Adjust the plot adding labels, legend, and saving in eps
    plt.xlabel("Temperature [K]")
    plt.ylabel(r"Volume [$\AA^3$]")
    plt.legend()
    plt.tight_layout()
    plt.savefig("thermal_expansion.png")
    plt.show()
    return 0

class  SrTiO3_free_energy_ab_initio(object):
    def __init__(self):
        # Initialize the DFT (Quantum Espresso) calculator for SrTiO3
        # The input data is a dictionary that encodes the pw.x input file namelist
        input_data = {
            'control' : {
                # Avoid writing wavefunctions on the disk
                'disk_io' : 'None',
                # Where to find the pseudopotential
                #'pseudo_dir' : '.'
            },
            'system' : {
                # Specify the basis set cutoffs
                'ecutwfc' : 60,   # Cutoff for wavefunction from 70
                'ecutrho' : 60*10, # Cutoff for the density ecutwfc*10
                # Information about smearing (it is a metal)
                'occupations' : 'fixed', # 'fixed' or 'smearing', smearing for conductors
                #'smearing' : 'mv',
                'degauss' : 0.03
            },
            'electrons' : {
                'conv_thr' : 1e-8,
                'conv_thr' : 1.e-09,
                'electron_maxstep' : 80,
                'mixing_beta' : 4.e-01
            }
        }

        # the pseudopotential for each chemical element
        # In this case O, Sr and Ti
        pseudopotentials = {'O' : 'O.pbesol-n-kjpaw_psl.1.0.0.UPF',
                            'Sr' : 'Sr.pbesol-spn-kjpaw_psl.1.0.0.UPF',
                            'Ti' : 'Ti.pbesol-spn-kjpaw_psl.1.0.0.UPF'}

        # the kpoints mesh and the offset
        kpts = (6,6,6)
        koffset = (0,0,0)

        # Specify the command to call quantum espresso
        command = 'mpirun -np 8 pw.x -i PREFIX.pwi > PREFIX.pwo'


        # Prepare the quantum espresso calculator
        self.calculator = CC.calculators.Espresso(input_data,
                                             pseudopotentials,
                                             command = command,
                                             kpts = kpts,
                                             koffset = koffset)

    def relax(self):

        TEMPERATURE = 300
        N_CONFIGS = 32
        MAX_ITERATIONS = 20
        START_DYN = 'harmonic_dyn'
        NQIRR = 10

        # Let us load the starting dynamical matrix
        SrTiO3_dyn = CC.Phonons.Phonons(START_DYN, NQIRR)
        SrTiO3_dyn.ForcePositiveDefinite()
        SrTiO3_dyn.Symmetrize()
        # Initialize the random ionic ensemble
        ensemble = sscha.Ensemble.Ensemble(SrTiO3_dyn, TEMPERATURE)

        # Initialize the free energy minimizer
        minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
        minim.set_minimization_step(0.01)

        # Initialize the NVT simulation
        mi_cluster = Send_to_cluster(hostname = 'diegom@ekhi.cfm.ehu.es',label = "SrTiO3_", n_pool = 5) #test with 5 pools for QE
        relax = sscha.Relax.SSCHA(minim, self.calculator, N_configs = N_CONFIGS,
                        max_pop = MAX_ITERATIONS, cluster = mi_cluster.cluster,
                        save_ensemble = True)

        # Define the I/O operations
        # To save info about the free energy minimization after each step
        ioinfo = sscha.Utilities.IOInfo()
        ioinfo.SetupSaving("minim_info")
        relax.setup_custom_functions(custom_function_post = ioinfo.CFP_SaveAll)


        # Run the NVT simulation (save the stress to compute the pressure)
        relax.relax(get_stress = True, sobol = True)

        # If instead you want to run a NPT simulation, use
        # The target pressure is given in GPa.
        #relax.vc_relax(target_press = 0)

        # You can also run a mixed simulation (NVT) but with variable lattice parameters
        #relax.vc_relax(fix_volume = True)

        # Now we can save the final dynamical matrix
        # And print in stdout the info about the minimization
        relax.minim.finalize()
        relax.minim.dyn.save_qe("sscha_T{}_dyn".format(TEMPERATURE))

def plot_dispersion_SrTiO3(PATH = "GX"):
    NQIRR = 10
    #CMAP = "Spectral_r"
    #PATH = "GX"
    #PATH = "GM"
    #PATH = "GR"

    N_POINTS = 1000

    SPECIAL_POINTS = {"G": [0,0,0],
                      "X": [0, 0, .5],
                      "M": [0, .5, .5],
                      "R": [.5, .5, .5]}

    # Load the harmonic and sscha phonons
    harmonic_dyn = CC.Phonons.Phonons('harmonic_dyn', NQIRR)
    sscha_dyn = CC.Phonons.Phonons('sscha_T300_dyn', NQIRR)

    # Get the band path
    qpath, data = CC.Methods.get_bandpath(harmonic_dyn.structure.unit_cell,
                                          PATH,
                                          SPECIAL_POINTS,
                                          N_POINTS)
    xaxis, xticks, xlabels = data # Info to plot correclty the x axis

    # Get the phonon dispersion along the path
    harmonic_dispersion = CC.ForceTensor.get_phonons_in_qpath(harmonic_dyn, qpath)
    sscha_dispersion = CC.ForceTensor.get_phonons_in_qpath(sscha_dyn, qpath)

    nmodes = harmonic_dyn.structure.N_atoms * 3

    # Plot the two dispersions
    plt.figure(dpi = 150)
    ax = plt.gca()

    for i in range(nmodes):
        lbl=None
        lblsscha = None
        if i == 0:
            lbl = 'Harmonic'
            lblsscha = 'SSCHA'

        ax.plot(xaxis, harmonic_dispersion[:,i], color = 'k', ls = 'dashed', label = lbl)
        ax.plot(xaxis, sscha_dispersion[:,i], color = 'r', label = lblsscha)

    # Plot vertical lines for each high symmetry points
    for x in xticks:
        ax.axvline(x, 0, 1, color = "k", lw = 0.4)
    ax.axhline(0, 0, 1, color = 'k', ls = ':', lw = 0.4)

    ax.legend()

    # Set the x labels to the high symmetry points
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

    ax.set_xlabel("Q path")
    ax.set_ylabel("Phonons [cm-1]")

    plt.tight_layout()
    plt.savefig("dispersion.png")
    plt.show()
    return 0

def main(args):
    Gold_calculation = Gold_free_energy()
    Gold_calculation.relax()
    Gold_calculation.harmonic()
    plot_dispersion()
    Gold_thermal_expansion = Thermal_expansion()
    Gold_thermal_expansion.temperature_cycle()
    plot_thermal_expansion()
    Gold_calculation_dft = Gold_free_energy_ab_initio()
    Gold_calculation_dft.relax()
    plot_dispersion()
    SrTiO3_calculation = SrTiO3_free_energy_ab_initio()
    SrTiO3_calculation.relax()
    plot_dispersion()
    #return 0
    #raise SystemExit
    #sys.exit()

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
