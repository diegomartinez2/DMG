#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  SSCHA_spectarl_ploter.py
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
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.ForceTensor

import numpy as np
import matplotlib.pyplot as plt

import sys,os

def plot_dispersion_SrTiO3(PATH ,sschaT300_file ,sschaT200_file ,sschaT100_file ,sschaT50_file ,sschaT0_file):
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
    sschaT300_dyn = CC.Phonons.Phonons(sschaT300_file, NQIRR)
    sschaT200_dyn = CC.Phonons.Phonons(sschaT200_file, NQIRR)
    sschaT100_dyn = CC.Phonons.Phonons(sschaT100_file, NQIRR)
    sschaT50_dyn = CC.Phonons.Phonons(sschaT50_file, NQIRR)
    sschaT0_dyn = CC.Phonons.Phonons(sschaT0_file, NQIRR)

    # Get the band path
    qpath, data = CC.Methods.get_bandpath(sschaT300_dyn.structure.unit_cell,
                                          PATH,
                                          SPECIAL_POINTS,
                                          N_POINTS)
    xaxis, xticks, xlabels = data # Info to plot correclty the x axis

    # Get the phonon dispersion along the path
    sschaT300_dispersion = CC.ForceTensor.get_phonons_in_qpath(sschaT300_dyn, qpath)
    sschaT200_dispersion = CC.ForceTensor.get_phonons_in_qpath(sschaT200_dyn, qpath)
    sschaT100_dispersion = CC.ForceTensor.get_phonons_in_qpath(sschaT100_dyn, qpath)
    sschaT50_dispersion = CC.ForceTensor.get_phonons_in_qpath(sschaT50_dyn, qpath)
    sschaT0_dispersion = CC.ForceTensor.get_phonons_in_qpath(sschaT0_dyn, qpath)

    nmodes = sschaT300_dyn.structure.N_atoms * 3

    # Plot the two dispersions
    plt.figure(dpi = 150)
    ax = plt.gca()

    for i in range(nmodes):
        lbl=None
        lblsscha = None
        if i == 0:
            #lbl = 'Harmonic'
            lblsschaT300 = 'SSCHA T=300K'
            lblsschaT200 = 'SSCHA T=200K'
            lblsschaT100 = 'SSCHA T=100K'
            lblsschaT50 = 'SSCHA T=50K'
            lblsschaT0 = 'SSCHA T=0K'
        else:
            lblsschaT300 = ''
            lblsschaT200 = ''
            lblsschaT100 = ''
            lblsschaT50 = ''
            lblsschaT0 = ''

        #ax.plot(xaxis, harmonic_dispersion[:,i], color = 'k', ls = 'dashed', label = lbl)
        ax.plot(xaxis, sschaT300_dispersion[:,i], color = 'r', label = lblsschaT300)
        ax.plot(xaxis, sschaT200_dispersion[:,i], color = 'y', label = lblsschaT200)
        ax.plot(xaxis, sschaT100_dispersion[:,i], color = 'k', label = lblsschaT100)
        ax.plot(xaxis, sschaT50_dispersion[:,i], color = 'b', label = lblsschaT50)
        ax.plot(xaxis, sschaT0_dispersion[:,i], color = 'r', linestyle = '--', label = lblsschaT0)

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

#    plt.legend()
    plt.tight_layout()
    plt.savefig("dispersion_all.png")
    plt.show()
    return 0

def main(args):
    plot_dispersion_SrTiO3(PATH = "GXMGRM",sschaT300_file = './T300/Hessian/dyn_hessian_', sschaT200_file='./200K/Hessian_200K/dyn_hessian_',sschaT100_file = './T100/Hessian/dyn_hessian_',sschaT50_file = './T50/Hessian/dyn_hessian_', sschaT0_file='./0K/Hessian_0K/dyn_hessian_')

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
