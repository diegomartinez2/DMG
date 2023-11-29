#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Spectral_plasmons_analysis.py
#
#  Copyright 2023 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
# ---------------------------
# Importación de los módulos
# ---------------------------
import Thermal
import Spectral_function_calc

from __future__ import print_function
from __future__ import division
# Import the modules to read the dynamical matrix
import cellconstructor as CC
import cellconstructor.Phonons
# Import the SCHA modules
import sscha, sscha.Ensemble

import numpy as np
import cellconstructor.ForceTensor
import cellconstructor.ThermalConductivity
import time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def main(args):
    """
    """
    if (len( sys.argv ) > 1):

    else:
        print ("Arguments are [...]. In that order, separated by simple spaces")
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
