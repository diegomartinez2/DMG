import cellconstructor as CC
import cellconstructor.ForceTensor
import cellconstructor.Structure
import cellconstructor.Spectral

import numpy as np

# Initialize the tensor3 object
# We need 2nd FCs of the used grid to configure the supercell.
# For example, we can use the sscha final auxiliary dynamical matrices
dyn = CC.Phonons.Phonons("dyn_end_population3_",3)
supercell = dyn.GetSupercell()
tensor3 = CC.ForceTensor.Tensor3(dyn.structure,
                                dyn.structure.generate_supercell(supercell),
                                supercell)

# Assign the tensor3 values
d3 = np.load("d3_realspace_sym.npy")*2.0 # The 2 factor is because of units, needs to be passed to Ry
tensor3.SetupFromTensor(d3)

# Center and apply ASR, which is needed to interpolate the third order force constant
tensor3.Center()
tensor3.Apply_ASR()

# Print the tensor if you want, uncommenting the next line
#tensor3.WriteOnFile(fname="FC3",file_format='D3Q')

# Calculate the spectral function at Gamma in the no-mode mixing approximation
# keeping the energy dependence on the self-energy.
#
# An interpolation grid needs to be used (and one needs to check convergence with
# respect to it considering different values of the smearing)


# interpolation grid
k_grid=[20,20,20]

#
G=[0.0,0.0,0.0]

CC.Spectral.get_diag_dynamic_correction_along_path(dyn=dyn,
                                                   tensor3=tensor3,
                                                   k_grid=k_grid,
                                                   q_path=G,
                                                   T = 200,                             # The temperature for the calculation
                                                   e1=145, de=0.1, e0=0,                # The energy grid in cm-1
                                                   sm1=1.0, nsm=1, sm0=1.0,             # The smearing \eta for the analytic continuation
                                                   filename_sp = 'nomm_spectral_func')  # Output file name

# Now perform the calculation of the spectral function in a
# path of q points where the list of q points is gicen in 2pi/a units, with
# a the lattice parameter given in Arnstrong

CC.Spectral.get_diag_dynamic_correction_along_path(dyn=dyn,
                                                   tensor3=tensor3,
                                                   k_grid=k_grid,
                                                   q_path_file="XGX.dat",
                                                   T = 200.0,
                                                   e1=145, de=0.1, e0=0,
                                                   sm1=1.0, nsm=1, sm0=1.0,
                                                   filename_sp = 'nomm_spectral_func_in_path')

def plot_dispersion_SrTiO3(PATH = "GXMGRX", NQIRR = 4, Temperatura = 300):

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
    harmonic_dyn = CC.Phonons.Phonons('harmonic_dyn', NQIRR)
    sscha_dyn = CC.Phonons.Phonons('sscha_T{}_dyn'.format(Temperatura), NQIRR)

    # Get the band path
    qpath, data = CC.Methods.get_bandpath(harmonic_dyn.structure.unit_cell,
                                          PATH,
                                          SPECIAL_POINTS,
                                          N_POINTS)
    xaxis, xticks, xlabels = data # Info to plot correclty the x axis

    # Get the phonon dispersion along the path
    harmonic_dispersion = CC.ForceTensor.get_phonons_in_qpath(harmonic_dyn, qpath)
    sscha_dispersion = CC.ForceTensor.get_phonons_in_qpath(sscha_dyn, qpath)

    nmodes = harmonic_dyn.structure.N_atoms * 3

    # Plot the two dispersions
    plt.figure(dpi = 150)
    ax = plt.gca()

    for i in range(nmodes):
        lbl=None
        lblsscha = None
        if i == 0:
            lbl = 'Harmonic'
            lblsscha = 'SSCHA'

        ax.plot(xaxis, harmonic_dispersion[:,i], color = 'k', ls = 'dashed', label = lbl)
        ax.plot(xaxis, sscha_dispersion[:,i], color = 'r', label = lblsscha)

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

    plt.tight_layout()
    plt.savefig("dispersion{}.png".format(PATH))
    plt.show()
    return 0
