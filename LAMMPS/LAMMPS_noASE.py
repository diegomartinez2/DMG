#!/usr/local/bin/python
import os
import lammps
import numpy as np
import pandas as pd
#Here we run simulation on a methane box using a Lenard Jones potential. Let’s write a LAMMPS input file in which you fill a box with 2000 methane particles using the TraPPE potential. The syntax is the same as the one of a classic LAMMPS input.
methane_inp = """# LAMMPS input file for Methane

units real
atom_style atomic

region myRegion block 0. 400. 0. 400. 0. 400. units box
create_box 1 myRegion
create_atoms 1 random 2000 03042019 myRegion

# LJ Methane: TraPPE potential
# eps / k_B = 148.0
# eps = 148.0 * k_B * N_a * 1e-3 / 4.184
#
pair_style lj/cut 14.0
pair_coeff   1   1    0.2941    3.730
pair_modify shift no mix arithmetic tail yes
mass 1 16.04

run_style verlet
neighbor 2.0 bin
neigh_modify delay 10

timestep 2.0
thermo_style multi
thermo 200
"""
#In order to run the calculations, you have to initialize the lammps object and execute LAMMPS commands through it.
# setup the lammps object and reduce verbosity
# lmp = lammps.lammps(cmdargs=["-log", "none", "-nocite"])
lmp = lammps.lammps(cmdargs=["-nocite"])
# load the initialization commands
lmp.commands_string(methane_inp)
#The simulations is still alive and you can for exemple, run one NVE step just to compute energy.
lmp.commands_list([
    "fix NVE all nve",
    "run 0"
])

etotal_0 = lmp.get_thermo("etotal")
print(f"E_total = {etotal_0:.2f} kcal.mol-1")
#Let’s now run a minimization and get back energies:
lmp.commands_string("minimize 1.0e-4 1.0e-6 1000 1000")
etotal_min = lmp.get_thermo("etotal")

print(f"E_total after minimization {etotal_min:.2f} kcal.mol-1")
print(f"Delta E_total = {etotal_min - etotal_0:.2f} kcal.mol-1")
#Let’s now set up a short NVT simulations.
lmp.commands_list([
    "velocity all create 300. 03042019 dist gaussian mom yes rot yes",
    "fix NVT all nvt temp 300. 300. $(100.0 * dt)",
    "dump trj all custom 10 traj.lammpstrj id type element x y z",
    "run 1000",
])
#If you want to deal with the coordinates of the methane molecules you can extract them in a numpy array:
coords = lmp.numpy.extract_atom("x")
print(coords.shape)
#At the end, close the LAMMPS execution:
lmp.close()
