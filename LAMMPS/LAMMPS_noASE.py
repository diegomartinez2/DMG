import os
import lammps
import numpy as np
import pandas as pd

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

# setup the lammps object and reduce verbosity
# lmp = lammps.lammps(cmdargs=["-log", "none", "-nocite"])
lmp = lammps.lammps(cmdargs=["-nocite"])
# load the initialization commands
lmp.commands_string(methane_inp)
