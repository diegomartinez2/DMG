"""
***Environment variables required in your shell:***

#!/usr/bin/bash

# ASE LAMMPSRun
export LAMMPS_COMMAND="/usr/bin/mpirun -np 4 lmp_mpi"

# LAMMPS
export LAMMPSPATH="$HOME/lib/lammps-ro"
export PYTHONPATH="$LAMMPSPATH/python:$PYTHONPATH"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$LAMMPSPATH/src"

***Installation script:***

#!/usr/bin/bash
git clone git://git.lammps.org/lammps-ro.git $LAMMPSPATH
cd $LAMMPSPATH/src
make yes-python
make yes-manybody
make mpi
make mpi mode=shlib

"""

#!/usr/local/bin/python
"""
LAMMPS calculation with ASE for C cristal (diamond)
"""
from ase import Atoms, Atom
from ase.calculators.lammpsrun import LAMMPS
from ase.visualize import view
import numpy as np
import os

# LAMMPS context information
lmp_path = os.getenv("LAMMPSPATH")
potential = os.path.join(lmp_path, "potentials/CH.airebo")
files = [potential]
parameters = {"mass": ["* 1.0"],
                "pair_style": "airebo 6.5 1 1",
                "pair_coeff": ['* * ' + potential + ' C']}
calc = LAMMPS(parameters=parameters, files=files)

# Graphene structure
a = 2.46
a1 = a * np.array([3.0**0.5/2., -1./2., 0.])
a2 = a * np.array([3.0**0.5/2., 1./2., 0.])
a3 = np.array([0., 0., 18.])
atoms = Atoms([Atom('C', 1./2. * a3),
                Atom('C', 1./3. * a1 + 1./3. * a2 + 1./2. * a3)],
                cell=[a1, a2, a3], pbc=True)

# Calculate the energy
atoms.set_calculator(calc)
energy = atoms.get_potential_energy()

print("Energy: {:0.3f} eV".format(energy))
