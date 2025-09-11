#!/bin/sh
./lmp_mpi_chimes -in lammps_relax.in
python VASP_LAMMPS_converter2.py
