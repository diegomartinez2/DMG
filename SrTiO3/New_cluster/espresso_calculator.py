import cellconstructor as CC
from cellconstructor.calculators import Espresso

import sys, os
import numpy as np

def get_calculator():

    pseudopotentials = {
    'O' : 'O.pbesol-n-kjpaw_psl.1.0.0.UPF',
    'Sr' : 'Sr.pbesol-spn-kjpaw_psl.1.0.0.UPF',
    'Ti' : 'Ti.pbesol-spn-kjpaw_psl.1.0.0.UPF'
    #"H" : "H.pbe-rrkjus_psl.1.0.0.UPF"
    }

    imput_params = {
        "control" : {
            "tprnfor" : True,
            "tstress" : True,
            "disk_io" : "none" # Avoid writing wavefunctions on the disk
            },
        "system" : {
            "degauss" : 0.03,
            #"smearing" : "mv",# remove for non-conductors and put occupations to 'fixed'
            "ecutwfc" : 35,# Cutoff for wavefunction
            "ecutrho" : 350,# Cutoff for the density
            "occupations" : "fixed"# 'fixed' or 'smearing', smearing for conductors
            },
        "electrons" : {
            #"conv_thr" : 1e-10
            #'conv_thr' : 1e-8,
            'conv_thr' : 1.e-09,
            'electron_maxstep' : 80,
            'mixing_beta' : 4.e-01
            }
    }

    kpts = (4,4,4)
    koffset = (1, 1, 1)

    calc = Espresso(pseudopotentials = pseudopotentials,
                    input_data = input_params,
                    kpts = kpts,
                    koffset = koffset)

    return calc
