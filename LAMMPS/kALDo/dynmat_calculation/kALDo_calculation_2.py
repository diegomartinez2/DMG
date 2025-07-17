#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# ASE imports
from ase.build import bulk, surface
from ase.build import cut
from ase.io import write, lammpsdata
from ase.visualize.plot import plot_atoms # Specific import for plotting atoms

# kALDo imports
from kaldo.forceconstants import ForceConstants
from kaldo.phonons import Phonons
from kaldo.conductivity import Conductivity
import kaldo.controllers.plotter as plotter # Alias for convenience

# --- Configuration ---
# Geometry and Supercell
NX, NY, NZ = 2, 2, 2
print(f"Supercell dimensions: {NX}x{NY}x{NZ}")

# k-point grid for phonon calculations
KX, KY, KZ = 10, 10, 10
KPTS = [KX, KY, KZ]

# Simulation Temperature
TEMPERATURE = 300 # Kelvin
IS_CLASSIC = False # Quantum or Classic statistics

# Material specific parameters (for 2D materials, otherwise Lz/dist = 1)
IS_2D_MATERIAL = False # Set to True for 2D materials
Lz = 1.0 # Box size in z-direction in LAMMPS (if applicable for 2D)
DIST = 1.0 # Thickness of the 2D material (if applicable)

# --- Plotting Configuration ---
# Denote plot default format (using plt.rcParams for consistency)
LINE_WIDTH_AXES = 2
FONT_SIZE = 12
TICK_LENGTH_MAJOR = 4
TICK_WIDTH_MAJOR = 2
TICK_LENGTH_MINOR = 2

plt.rcParams['font.size'] = FONT_SIZE
plt.rcParams['axes.linewidth'] = LINE_WIDTH_AXES
plt.rcParams['text.usetex'] = False # Set to True if LaTeX is installed and needed for math text
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['mathtext.fontset'] = 'cm' # or 'stixsans'

# --- Helper Functions ---

def cumulative_cond_cal(observables, kappa_tensor, prefactor=1/3):
    """
    Compute cumulative conductivity based on either frequency or mean-free path.

    Parameters
    ----------
    observables : np.ndarray
        Either phonon frequency or mean-free path.
    kappa_tensor : np.ndarray
        Conductivity tensor (e.g., from kALDo, shape (modes, 3, 3)).
    prefactor : float, optional
        Prefactor to average kappa tensor, 1/3 for bulk material.
        Default is 1/3. For 2D, this may need adjustment or be omitted.

    Returns
    -------
    observables_sorted : np.ndarray
        Sorted phonon frequency or mean-free path.
    cumulative_kappa : np.ndarray
        Cumulative conductivity.
    """
    # Sum over kappa by directions (e.g., diagonal components for isotropic materials)
    # kappa_tensor has shape (n_modes, 3, 3)
    # We want sum over diagonal components for average conductivity
    # np.einsum('maa->m', prefactor * kappa_tensor) sums over the second and third
    # dimensions for each mode, taking the trace of each 3x3 submatrix.
    # If prefactor is 1/3, this calculates (k_xx + k_yy + k_zz) / 3 for each mode.
    kappa = np.einsum('maa->m', prefactor * kappa_tensor)

    # Sort observables and cumulative sum kappa
    observables_argsort_indices = np.argsort(observables)
    cumulative_kappa = np.cumsum(kappa[observables_argsort_indices])
    return observables[observables_argsort_indices], cumulative_kappa

def set_plot_style(ax_list):
    """
    Applies consistent styling to a list of matplotlib axes.
    """
    for ax in ax_list:
        ax.tick_params(which='major', length=TICK_LENGTH_MAJOR, width=TICK_WIDTH_MAJOR)
        ax.tick_params(which='minor', length=TICK_LENGTH_MINOR, width=TICK_WIDTH_MAJOR)
        ax.tick_params(which='both', axis='both', direction='in',
                       right=True, top=True)

        # Set spine colors and linewidth (can be customized)
        ax.spines['bottom'].set_color('black')
        ax.spines['top'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.spines['right'].set_color('black')

        ax.spines['bottom'].set_linewidth(LINE_WIDTH_AXES)
        ax.spines['top'].set_linewidth(LINE_WIDTH_AXES)
        ax.spines['left'].set_linewidth(LINE_WIDTH_AXES)
        ax.spines['right'].set_linewidth(LINE_WIDTH_AXES)

        for t in ax.xaxis.get_ticklines():
            t.set_color('black')
            t.set_linewidth(LINE_WIDTH_AXES)
        for t in ax.yaxis.get_ticklines():
            t.set_color('black')
            t.set_linewidth(LINE_WIDTH_AXES)

def main():
    # --- 1. Clean previous results (optional, use with caution) ---
    # Uncomment if you want to remove old data. Be careful with 'rm -r *'
    # try:
    #     subprocess.run('rm -rf *fd* *xyz ge-* plots *lmp *log* graphene', shell=True, check=True)
    # except subprocess.CalledProcessError as e:
    #     print(f"Warning: Could not remove old files. Error: {e}")

    # --- 2. Load Atomic Structure ---
    print("\n--- Loading Atomic Structure ---")
    try:
        # Adjust Z_of_type according to your LAMMPS data file and atom types
        # Example: {1:16, 2:6, 3:29} means type 1 is S (16), type 2 is C (6), type 3 is Cu (29)
        atoms = lammpsdata.read_lammps_data('relax.dat', Z_of_type={1:16, 2:6, 3:29}, style='charge', units='metal')
        print("Successfully read LAMMPS relaxed data from 'relax.dat'.")
    except FileNotFoundError:
        print("Error: 'relax.dat' not found. Please ensure the relaxed LAMMPS data file exists.")
        return # Exit if essential file is missing

    # Create supercell and save to XYZ
    supercell_dims = np.array([NX, NY, NZ])
    replicated_atoms = atoms.copy().repeat(supercell_dims)
    write('replicated_atoms.xyz', format='extxyz', images=replicated_atoms)
    print(f"Saved replicated_atoms.xyz with supercell {NX}x{NY}x{NZ}.")

    # Plot initial atoms figure
    fig_atoms, ax_atoms = plt.subplots(figsize=(6, 6))
    plot_atoms(replicated_atoms, ax_atoms, radii=0.3, rotation=('0x,0y,0z'))
    ax_atoms.set_title(f'Replicated Atoms ({NX}x{NY}x{NZ} Supercell)')
    set_plot_style([ax_atoms]) # Apply consistent style
    fig_atoms.savefig("ase_atoms.png", dpi=300, bbox_inches='tight')
    plt.close(fig_atoms) # Close figure to free memory
    print("Saved 'ase_atoms.png'.")

    # --- 3. Set Force Constants with kALDo ---
    print("\n--- Setting Force Constants ---")
    # kALDo assumes FC files are in the specified format and folder.
    # 'Dyn.form' for 2nd order, 'THIRD' for 3rd order.
    # `only_second=True` means only 2nd order FCs are used.
    # If 3rd order is needed for QHGK or other calculations, remove `only_second=True`.
    try:
        forceconstants = ForceConstants.from_folder(
            folder='./', supercell=supercell_dims, only_second=True, format='lammps'
        )
        print("Force constants set from LAMMPS format data:")
        print(f"Folder: './', Config_filename: 'replicated_atoms.xyz', ASE_format: 'xyz'")
        print(f"2nd_order_FC: 'Dyn.form', 3rd_order_FC: 'THIRD' (not used with only_second=True)")
    except Exception as e:
        print(f"Error setting force constants: {e}")
        return

    # --- 4. Create Phonon Object ---
    print("\n--- Creating Phonon Object ---")
    k_label = f"{KX}_{KY}_{KZ}" # Consistent k-label
    phonon_folder = f'Supercell' # kALDo will create folders inside this

    phonons = Phonons(forceconstants=forceconstants,
                      kpts=KPTS,
                      is_classic=IS_CLASSIC,
                      temperature=TEMPERATURE,
                      folder=phonon_folder,
                      storage='numpy')
    print(f"Phonon object created for {k_label} k-grid at {TEMPERATURE}K.")

    # --- 5. Plot Phonon Dispersion and DOS ---
    print("\n--- Plotting Phonon Dispersion and DOS ---")

    # Define k-path for dispersion (adjust based on your material's symmetry)
    # Using the special points defined in ASE's bandpath for the specific cell.
    band_path = atoms.cell.bandpath(path='GXVRGY', npoints=200,
                                    special_points={'G': [0, 0, 0], 'X': [1/2, 0, 0], 'Y': [0, 1/2, 0],
                                                    'Z': [0, 0, 1/2], 'R': [1/2, 1/2, 1/2], 'T': [0, 1/2, 1/2],
                                                    'U': [1/2, 0, 1/2], 'V': [1/2, 1/2, 0]}) # V was [1/2, 0, 1/2] in original, corrected to [1/2,1/2,0] if meant for orthorhombic

    # Plot Dispersion
    fig_dispersion, ax_dispersion = plt.subplots(figsize=(8, 6))
    plotter.plot_dispersion(phonons, with_velocity=True, is_showing=False, # Don't show immediately
                            manually_defined_path=band_path, fig=fig_dispersion, ax=ax_dispersion)
    ax_dispersion.set_title('Phonon Dispersion Relation')
    set_plot_style([ax_dispersion])
    fig_dispersion.savefig("ase_dispersion.png", dpi=300, bbox_inches='tight')
    plt.close(fig_dispersion)
    print("Saved 'ase_dispersion.png'.")

    # Plot DOS
    fig_dos, ax_dos = plt.subplots(figsize=(8, 6))
    # `plotter.plot_dos` usually returns ax, figure directly, or plot to current ax.
    # If it creates a new figure, pass fig, ax arguments, or capture them.
    plotter.plot_dos(phonons, p_atoms=None, bandwidth=2, n_points=200, filename='dos',
                     is_showing=False, fig=fig_dos, ax=ax_dos)
    ax_dos.set_title('Phonon Density of States (DOS)')
    set_plot_style([ax_dos])
    fig_dos.savefig("ase_DOS.png", dpi=300, bbox_inches='tight')
    plt.close(fig_dos)
    print("Saved 'ase_DOS.png'.")

    # --- 6. Calculate Conductivity ---
    print("\n--- Calculating Conductivity (Inverse Method) ---")
    # Conductivity from inversion
    # The `conductivity` property returns a (modes, 3, 3) tensor by default.
    # Summing over axis=0 gives the total 3x3 conductivity matrix.
    total_cond_matrix = Conductivity(phonons=phonons, method='inverse').conductivity.sum(axis=0)

    # Adjust conductivity for 2D materials if applicable
    kappa_factor = (Lz / DIST) if IS_2D_MATERIAL and DIST != 0 else 1.0
    effective_cond_matrix = kappa_factor * total_cond_matrix

    # For 3D materials, mean of diagonal elements (xx, yy, zz)
    mean_conductivity_3D = np.mean(np.diag(effective_cond_matrix))
    print(f'Total Conductivity (W/m-K): {mean_conductivity_3D:.3f}')
    print("Full conductivity matrix:")
    print(effective_cond_matrix)

    # --- 7. Load and Process Phonon Data for Detailed Plots ---
    print("\n--- Loading Phonon Data for Detailed Plots ---")
    data_folder = os.path.join(phonon_folder, k_label)

    # Check if data files exist
    required_files = ['frequency.npy', f'{TEMPERATURE}/quantum/heat_capacity.npy', 'velocity.npy',
                      f'{TEMPERATURE}/quantum/_ps_and_gamma.npy', f'{TEMPERATURE}/quantum/bandwidth.npy']
    for f_name in required_files:
        if not os.path.exists(os.path.join(data_folder, f_name)):
            print(f"Error: Missing data file: {os.path.join(data_folder, f_name)}. Please check if kALDo calculation completed successfully.")
            return

    frequency = np.load(os.path.join(data_folder, 'frequency.npy'), allow_pickle=True)
    cv = np.load(os.path.join(data_folder, f'{TEMPERATURE}/quantum/heat_capacity.npy'), allow_pickle=True)
    group_velocity = np.load(os.path.join(data_folder, 'velocity.npy'))
    #[:,0] to select the phase space (first column)
    phase_space = np.load(os.path.join(data_folder, f'{TEMPERATURE}/quantum/_ps_and_gamma.npy'), allow_pickle=True)[:, 0]

    # Compute norm of group velocity
    # Convert unit from Angstrom/Picosecond to km/s (1 A/ps = 0.1 km/s)
    group_velocity_norm = np.linalg.norm(group_velocity.reshape(-1, 3), axis=1) / 10.0
    print("Group velocity norm computed.")

    # --- 8. Plot Observables (Cv, Group Velocity, Phase Space) ---
    print("--- Plotting Phonon Observables ---")
    fig1, axes1 = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    set_plot_style(axes1) # Apply style to all subplots

    # Plot Cv vs Frequency
    axes1[0].scatter(frequency.flatten(), 1e23 * cv.flatten(),
                     facecolor='w', edgecolor='r', s=10, marker='8')
    axes1[0].set_ylabel(r"$C_{v}$ ($10^{23}$ J/K)", fontsize=FONT_SIZE)
    axes1[0].set_xlabel('Frequency (THz)', fontsize=FONT_SIZE)
    # Set y-limits carefully to avoid empty plot if data is sparse or small
    if cv.flatten()[3:].size > 0: # Check if there's data to plot from index 3
        cv_min = 1e23 * cv.flatten()[3:].min()
        cv_max = 1e23 * cv.flatten()[3:].max()
        if cv_max > cv_min: # Avoid division by zero or invalid range
            axes1[0].set_ylim(0.9 * cv_min, 1.05 * cv_max)
        else: # Handle cases where all values might be identical
            axes1[0].set_ylim(cv_min * 0.9, cv_min * 1.1 + 1e-10) # Small offset
    else: # If no data points, set a default small range or warning
        axes1[0].set_ylim(0, 1e23 * 1e-3) # Example default range if no data

    # Plot Group Velocity vs Frequency
    axes1[1].scatter(frequency.flatten(), group_velocity_norm,
                     facecolor='w', edgecolor='r', s=10, marker='^')
    axes1[1].set_xlabel('Frequency (THz)', fontsize=FONT_SIZE)
    axes1[1].set_ylabel(r'$|v| \ (\rm{km/s})$', fontsize=FONT_SIZE) # Use \rm for roman text

    # Plot Phase Space vs Frequency
    axes1[2].scatter(frequency.flatten(), phase_space,
                     facecolor='w', edgecolor='r', s=10, marker='o')
    axes1[2].set_xlabel('Frequency (THz)', fontsize=FONT_SIZE)
    axes1[2].set_ylabel('Phase Space', fontsize=FONT_SIZE)

    fig1.savefig("ase_multiplot.png", dpi=300, bbox_inches='tight')
    plt.close(fig1)
    print("Saved 'ase_multiplot.png'.")

    # --- 9. Load and Process Scattering Rate and Mean Free Path ---
    print("\n--- Loading Scattering Rate and Mean Free Path ---")
    scattering_rate = np.load(os.path.join(data_folder, f'{TEMPERATURE}/quantum/bandwidth.npy'), allow_pickle=True)

    # Check if scattering_rate is an array or a single value
    if scattering_rate.ndim > 0 and scattering_rate.size > 0:
        # Avoid division by zero for lifetime calculation
        life_time = np.zeros_like(scattering_rate, dtype=float)
        non_zero_rates = scattering_rate != 0
        life_time[non_zero_rates] = 1.0 / scattering_rate[non_zero_rates]
    else: # Handle scalar or empty array cases
        life_time = np.array([0.0]) # Or handle as error
        if scattering_rate == 0:
            print("Warning: Scattering rate is zero. Lifetime will be infinite/undefined.")
        else:
            life_time = 1.0 / scattering_rate.astype(float)

    print(f"Scattering rate (first few values): {scattering_rate.flatten()[:5]}")
    print(f"Lifetime (first few values): {life_time.flatten()[:5]}")

    # Load mean free path for each direction
    mean_free_path_list = []
    for i in range(3):
        mfp_file = os.path.join(data_folder, f'{TEMPERATURE}/quantum/inverse/mean_free_path_{i}.dat')
        if not os.path.exists(mfp_file):
            print(f"Error: Missing mean free path file: {mfp_file}")
            return
        mean_free_path_list.append(np.loadtxt(mfp_file))

    mean_free_path = np.array(mean_free_path_list).T
    # Convert unit from Angstrom to Nanometer (1 A = 0.1 nm)
    mean_free_path_norm = np.linalg.norm(mean_free_path.reshape(-1, 3), axis=1) / 10.0
    print("Mean free path norm computed.")

    # --- 10. Plot Lifetime, Scattering Rate, Mean Free Path ---
    print("--- Plotting Lifetime, Scattering Rate, Mean Free Path ---")
    fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    set_plot_style(axes2)

    # Plot Lifetime vs Frequency
    axes2[0].scatter(frequency.flatten(), life_time,
                     facecolor='w', edgecolor='r', s=10, marker='s')
    axes2[0].set_yscale('log')
    axes2[0].set_xscale('log')
    axes2[0].set_ylabel(r'$\tau \ (\rm{ps})$', fontsize=FONT_SIZE)
    axes2[0].set_xlabel('Frequency (THz)', fontsize=FONT_SIZE)

    # Plot Scattering Rate vs Frequency
    axes2[1].scatter(frequency.flatten(), scattering_rate,
                     facecolor='w', edgecolor='r', s=10, marker='d')
    axes2[1].set_ylabel(r'$\Gamma \ (\rm{THz})$', fontsize=FONT_SIZE)
    axes2[1].set_xlabel('Frequency (THz)', fontsize=FONT_SIZE)

    # Plot Mean Free Path vs Frequency
    axes2[2].scatter(frequency.flatten(), mean_free_path_norm,
                     facecolor='w', edgecolor='r', s=10, marker='8')
    axes2[2].set_ylabel(r'$\lambda \ (\rm{nm})$', fontsize=FONT_SIZE)
    axes2[2].set_xlabel('Frequency (THz)', fontsize=FONT_SIZE)
    axes2[2].set_yscale('log')
    axes2[2].set_ylim([1e-1, 1e5]) # Ensure reasonable limits

    fig2.savefig("ase_multiplot2.png", dpi=300, bbox_inches='tight')
    plt.close(fig2)
    print("Saved 'ase_multiplot2.png'.")

    # --- 11. Plot Cumulative Conductivity ---
    print("\n--- Plotting Cumulative Conductivity ---")
    kappa_tensor_modes_3x3 = np.zeros((mean_free_path.shape[0], 3, 3))
    for i in range(3):
        for j in range(3):
            cond_file = os.path.join(data_folder, f'{TEMPERATURE}/quantum/inverse/conductivity_{i}_{j}.dat')
            if not os.path.exists(cond_file):
                print(f"Error: Missing conductivity file: {cond_file}")
                return
            kappa_tensor_modes_3x3[:, i, j] = np.loadtxt(cond_file)

    # Apply factor for 2D materials (if applicable, otherwise factor is 1.0)
    kappa_tensor_modes_3x3_scaled = kappa_factor * kappa_tensor_modes_3x3

    # Sum over the 0th dimension to recover 3-by-3 total kappa matrix
    kappa_matrix_total = kappa_tensor_modes_3x3_scaled.sum(axis=0)

    print(f"Bulk thermal conductivity (W m^-1 K^-1):")
    print(kappa_matrix_total)

    # For 2D materials, it's common to average only in-plane components (xx, yy)
    if IS_2D_MATERIAL:
        mean_in_plane_kappa = np.mean(np.diag(kappa_matrix_total[0:2, 0:2]))
        print(f"Mean in-plane thermal conductivity (W m^-1 K^-1): {mean_in_plane_kappa:.1f}")
    else:
        mean_total_kappa = np.mean(np.diag(kappa_matrix_total))
        print(f"Mean total thermal conductivity (W m^-1 K^-1): {mean_total_kappa:.1f}")

    # Compute kappa per mode and cumulative representations
    # kappa_per_mode sums the diagonal contributions for each mode: k_xx + k_yy + k_zz
    kappa_per_mode_sum_diag = np.einsum('mdd->m', kappa_tensor_modes_3x3_scaled)

    # Use the helper function for cumulative conductivity
    freq_sorted, kappa_cum_wrt_freq = cumulative_cond_cal(
        frequency.flatten(), kappa_tensor_modes_3x3_scaled, prefactor=1/3 if not IS_2D_MATERIAL else 1.0/3.0 # Factor for average
    )
    lambda_sorted, kappa_cum_wrt_lambda = cumulative_cond_cal(
        mean_free_path_norm, kappa_tensor_modes_3x3_scaled, prefactor=1/3 if not IS_2D_MATERIAL else 1.0/3.0
    )

    fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    set_plot_style(axes3)

    # Plot kappa per mode vs Frequency
    axes3[0].scatter(frequency.flatten(), kappa_per_mode_sum_diag,
                     facecolor='w', edgecolor='r', s=10, marker='>')
    axes3[0].axhline(y=0, color='k', ls='--', lw=1)
    axes3[0].set_ylabel(r'$\kappa_{\rm{per\,mode}}\;\left(\frac{\rm{W}}{\rm{m}\cdot\rm{K}}\right)$', fontsize=FONT_SIZE)
    axes3[0].set_xlabel('Frequency (THz)', fontsize=FONT_SIZE)
    # axes3[0].set_ylim([-0.5, 6]) # Uncomment and adjust if specific limits are needed

    # Plot cumulative kappa vs Frequency
    axes3[1].plot(freq_sorted, kappa_cum_wrt_freq, 'r-')
    axes3[1].set_ylabel(r'$\kappa_{\rm{cumulative}, \omega}\;\left(\frac{\rm{W}}{\rm{m}\cdot\rm{K}}\right)$', fontsize=FONT_SIZE)
    axes3[1].set_xlabel('Frequency (THz)', fontsize=FONT_SIZE)
    axes3[1].legend(['Cumulative Kappa'], loc='lower right', frameon=False, fontsize=FONT_SIZE * 0.8)

    # Plot cumulative kappa vs Mean Free Path
    axes3[2].plot(lambda_sorted, kappa_cum_wrt_lambda, 'r-')
    axes3[2].set_xlabel(r'$\lambda \ (\rm{nm})$', fontsize=FONT_SIZE)
    axes3[2].set_ylabel(r'$\kappa_{\rm{cumulative}, \lambda}\;\left(\frac{\rm{W}}{\rm{m}\cdot\rm{K}}\right)$', fontsize=FONT_SIZE)
    axes3[2].set_xscale('log')
    axes3[2].set_xlim(1e-1,)

    fig3.savefig("ase_multiplot3.png", dpi=300, bbox_inches='tight')
    plt.close(fig3)
    print("Saved 'ase_multiplot3.png'.")

    print("\n--- End of code ---")

if __name__ == '__main__':
    main()
