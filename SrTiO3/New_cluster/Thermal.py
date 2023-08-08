from __future__ import print_function
from __future__ import division

import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.ForceTensor
import cellconstructor.ThermalConductivity
import time

dyn_prefix = 'final_dyn'
nqirr = 8

SSCHA_TO_MS = cellconstructor.ThermalConductivity.SSCHA_TO_MS
RY_TO_THZ = cellconstructor.ThermalConductivity.SSCHA_TO_THZ
dyn = CC.Phonons.Phonons(dyn_prefix, nqirr)

supercell = dyn.GetSupercell()

fc3 = CC.ForceTensor.Tensor3(dyn.structure,
dyn.structure.generate_supercell(supercell), supercell)

d3 = np.load("d3.npy")
fc3.SetupFromTensor(d3)
fc3 = CC.ThermalConductivity.centering_fc3(fc3)

mesh = [10,10,10]
smear = 0.03/RY_TO_THZ

tc = CC.ThermalConductivity.ThermalConductivity(dyn, fc3,
kpoint_grid = mesh, scattering_grid = mesh, smearing_scale = None,
smearing_type = 'constant', cp_mode = 'quantum', off_diag = False)

temperatures = np.linspace(200,1200,10,dtype=float)
start_time = time.time()
tc.setup_harmonic_properties(smear)
tc.write_harmonic_properties_to_file()

tc.calculate_kappa(mode = 'SRTA', temperatures = temperatures,
write_lifetimes = True, gauss_smearing = True, offdiag_mode = 'wigner',
kappa_filename = 'Thermal_conductivity_SRTA', lf_method = 'fortran-LA')

print('Calculated SSCHA kappa in: ', time.time() - start_time)

tc.calculate_kappa(mode = 'GK', write_lineshapes=False,
ne = 1000, temperatures = temperatures,
kappa_filename = 'Thermal_conductivity_GK')

print('Calculated SSCHA kappa in: ', time.time() - start_time)
# Save ThermalConductivity object for postprocessing.
tc.save_pickle()

#----------------------------------------------------------------
#from __future__ import print_function
#from __future__ import division

# Import the modules to read the dynamical matrix
#import numpy as np
#import cellconstructor as CC
#import cellconstructor.Phonons
#import cellconstructor.ForceTensor
#import cellconstructor.ThermalConductivity
#import time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

tc = CC.ThermalConductivity.load_thermal_conductivity()

# See at which temperatures we calculated stuff
tc.what_temperatures()

key = list(tc.lineshapes.keys()) # Get Ts for lineshapes

# DOS calculated from auxiliary force constants
harm_dos = tc.get_dos()
# Temperature dependent DOS calculated from lineshapes
# first two arrays are raw data
# second two is gaussian smoothed results \
#for the distance between energy points de
anharm_dos = tc.get_dos_from_lineshapes(float(key[-1]), de = 0.1)

def plot(self, harm_dos, anharm_dos):
    """# Plot results"""
    fig = plt.figure(figsize=(6.4, 4.8))
    gs1 = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs1[0, 0])
    ax.plot(harm_dos[0], harm_dos[1], 'k-',
    zorder=0, label = 'Harmonic')
    ax.plot(anharm_dos[0], anharm_dos[1], 'r-',
    zorder=0, label = 'Anharmonic raw @ ' + key[-1] + ' K')
    ax.plot(anharm_dos[2], anharm_dos[3], 'b-',
    zorder=0, label = 'Anharmonic smooth @ ' + key[-1] + ' K')
    ax.set_xlabel('Frequency (THz)')
    ax.set_ylabel('Density of states')
    ax.legend(loc = 'upper right')
    ax.set_ylim(bottom = 0.0)
    fig.savefig('test.pdf')
    plt.show()
    pass

def main(args):
    pass

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
