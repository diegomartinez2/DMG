#!/bin/bash
python extract.py --LAMMPS=relax.dat XFSET.harm* > DFSET_harmonic
python extract.py --LAMMPS=relax.dat XFSET.cubic* > DFSET_cubic
cat DFSET_harmonic DFSET_cubic > DFSET_merged
