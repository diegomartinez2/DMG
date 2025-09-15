#!/bin/sh
python displace.py --LAMMPS=LAMMPS-relax.dat --mag=0.01 --prefix harm  -pf My_displacements_patterns.pattern_HARMONIC
python displace.py --LAMMPS=LAMMPS-relax.dat --mag=0.04 --prefix cubic  -pf My_displacements_patterns.pattern_ANHARM3
