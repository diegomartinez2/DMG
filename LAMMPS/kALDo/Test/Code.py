#!/usr/local/bin/python
"""
kALDo calculation script, requires LAMMPS output file with...(dynamicak matrix or lammps output file???)
requires ASE, pylab, subprocess, and kALDo
"""
from ase.build import bulk
from ase.io import lammpsdata
from ase.calculators.lammpslib import LAMMPSlib
from kaldo.forceconstants import ForceConstants
from kaldo.conductivity import Conductivity
from kaldo.phonons import Phonons
from pylab import *
from ase.io import read,write
import subprocess

import kaldo.controllers.plotter as plotter
import matplotlib
import numpy as np
import warnings
warnings.filterwarnings("ignore")

#Functions
def cumulative_cond_cal(freq,full_cond,n_phonons):

  conductivity = np.einsum('maa->m', 1/3 * full_cond)
  conductivity = conductivity.reshape(n_phonons)
  cumulative_cond = np.zeros_like(conductivity)
  freq_reshaped = freq.reshape(n_phonons)

  for mu in range(cumulative_cond.size):
      single_cumulative_cond = conductivity[(freq_reshaped < freq_reshaped[mu])].sum()
      cumulative_cond[mu] = single_cumulative_cond

  return cumulative_cond

def cumulative_cond_cal(observables, kappa_tensor, prefactor=1/3):

    """Compute cumulative conductivity based on either frequnecy or mean-free path

       input:
       observables: (ndarray) either phonon frequnecy or mean-free path
       cond_tensor: (ndarray) conductivity tensor
       prefactor: (float) prefactor to average kappa tensor, 1/3 for bulk material

       ouput:
       observeables: (ndarray) sorted phonon frequency or mean-free path
       kappa_cond (ndarray) cumulative conductivity

    """

    # Sum over kappa by directions
    kappa = np.einsum('maa->m', prefactor * kappa_tensor)

    # Sort observables
    observables_argsort_indices = np.argsort(observables)
    cumulative_kappa = np.cumsum(kappa[observables_argsort_indices])
    return observables[observables_argsort_indices], cumulative_kappa

def set_fig_properties(ax_list, panel_color_str='black', line_width=2):
    tl = 4
    tw = 2
    tlm = 2

    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='in',
                       right=True, top=True)
        ax.spines['bottom'].set_color(panel_color_str)
        ax.spines['top'].set_color(panel_color_str)
        ax.spines['left'].set_color(panel_color_str)
        ax.spines['right'].set_color(panel_color_str)

        ax.spines['bottom'].set_linewidth(line_width)
        ax.spines['top'].set_linewidth(line_width)
        ax.spines['left'].set_linewidth(line_width)
        ax.spines['right'].set_linewidth(line_width)

        for t in ax.xaxis.get_ticklines(): t.set_color(panel_color_str)
        for t in ax.yaxis.get_ticklines(): t.set_color(panel_color_str)
        for t in ax.xaxis.get_ticklines(): t.set_linewidth(line_width)
        for t in ax.yaxis.get_ticklines(): t.set_linewidth(line_width)

#main
def main(args):
    """
    Test code for a kALDo calculation
    """
    atom_type={}

    for i in range(1,17):
        atom_type[i] = 6
    for i in range(17,27):
        atom_type[i] = 1
    for i in range(27,31):
        atom_type[i] = 16

    atoms=lammpsdata.read_lammps_data('T4_relaxed.lammps',Z_of_type=atom_type,sort_by_id='True',units = 'real',style = 'full')

    b1 = atoms.get_reciprocal_cell()[:,0]
    b2 = atoms.get_reciprocal_cell()[:,1]
    b3 = atoms.get_reciprocal_cell()[:,2]

    b1_length = np.linalg.norm(b1)
    b2_length = np.linalg.norm(b2)
    b3_length = np.linalg.norm(b3)
    # Calculate the number of divisions along each direction
    # N1 = 3
    # N2 = 4
    # N3 = 2

    # res1 = b1_length/N1
    # res2 = b2_length/N2
    # res3 = b3_length/N3

    # print(res1,res2,res3)

    res = 0.03
    N1 = b1_length/res
    N2 = b2_length/res
    N3 = b3_length/res

    print(np.round(N1),np.round(N2),np.round(N3))

    nrep = 3
    supercell = np.array([nrep, nrep, 1])

    write('replicated_atoms.xyz',format ='extxyz',images=atoms.copy().repeat(supercell))

    # Print reminder information
    print('Supercell structures and LAMMPS input generated.')
    print('Supercell dimension is: ' + str(supercell))

    forceconstants = ForceConstants.from_folder(folder='./',supercell=supercell,format='lammps')

    kx = int(np.round(N1))
    ky = int(np.round(N2))
    kz = int(np.round(N3))
    kpts = [kx, ky, kz]
    temperature = 10
    is_classic = False
    k_label = str(kx) + '_' + str(ky) + '_' + str(kz)
    folder='dynmat-' + k_label
    # Create a phonon object
    phonons = Phonons(forceconstants=forceconstants,
                    kpts=kpts,
                    is_classic=is_classic,
                    temperature=temperature,
                    folder=folder,
                    # is_unfolding=True,
                    storage='numpy')
    plotter.plot_dos(phonons,bandwidth=0.25,n_points=500)

    # Plot dispersion relation and group velocity in each direction
    kpath=atoms.cell.bandpath(path='GYCEGZDAEZ',npoints=500,special_points={'G': [0,0,0], 'Y': [0,0.5,0], 'Z': [0,0,0.5], 'A': [-0.5,0.5,0], 'E': [-0.5,0,0], 'C': [-0.5,0,0.5], 'D': [-0.5,0.5,0.5], 'H': [0,0.5,0.5]})
    plotter.plot_dispersion(phonons,manually_defined_path=kpath)
    # plotter.plot_dispersion(phonons,n_k_points=200)
    print('\n')

    # Denote data path
    data_folder = "./"

    # Derive symbols for high symmetry directions in Brillouin zone
    dispersion = np.loadtxt(data_folder  + 'plots/'+k_label+'/dispersion')
    point_names = np.loadtxt(data_folder  + 'plots/'+k_label+'/point_names', dtype=str)

    point_names_list = []
    for point_name in point_names:
        if point_name == 'G':
            point_name = r'$\Gamma$'
        elif point_name == 'U':
            point_name = 'U=K'
        point_names_list.append(point_name)

    # Load in the "tick-mark" values for these symbols
    q = np.loadtxt(data_folder  + 'plots/'+k_label+'/q')
    Q = np.loadtxt(data_folder  + 'plots/'+k_label+'/Q_val')
    dos_Si = np.load(data_folder + 'plots/'+k_label+'/dos.npy')

    fig = plt.figure(figsize=(4,3))
    set_fig_properties([gca()])
    plot(q[0], dispersion[0, 0], 'r-', ms=1)
    plot(q, dispersion, 'r-', ms=1)
    for i in range(1, 4):
        axvline(x=Q[i], ymin=0, ymax=2, ls='--',  lw=2, c="k")
    ylabel('Frequency (THz)', fontsize=14)
    xlabel(r'Wave vector ($\frac{2\pi}{a}$)', fontsize=14)
    gca().set_yticks(np.arange(0, 30, 5))
    xticks(Q, point_names_list)
    ylim([0, 18])
    xlim([Q[0], Q[-1]])
    dosax = fig.add_axes([0.91, .11, .17, .77])
    set_fig_properties([gca()])
    for d in np.expand_dims(dos_Si[1],0):
        dosax.plot(d, dos_Si[0],c='r')
    dosax.set_yticks([])
    dosax.set_xticks([])
    dosax.set_xlabel("DOS")
    ylim([0, 18])
    show()

    frequency = phonons.frequency.flatten(order='C')
    heat_capacity = phonons.heat_capacity.flatten(order='C')
    bandwidth = phonons.bandwidth.flatten(order='C')

    plt.figure()
    plt.scatter(frequency[3:], 1e23 * heat_capacity[3:], s=5)  # Get rid of the first three non-physical modes while plotting
    plt.xlabel("$\\nu$ (THz)", fontsize=16)
    plt.ylabel(r"$C_{v} \ (10^{23} \ J/K)$", fontsize=16)

    # Plot phonon bandwidth vs frequency
    print('\n')
    plt.plot(frequency[3:],bandwidth[3:],'.',markersize=10,label= 'broadening shape: ' + str(phonons.broadening_shape))
    plt.ylabel('$\Gamma$ (THz)', fontsize=25, fontweight='bold')
    plt.xlabel("$\\nu$ (THz)", fontsize=25, fontweight='bold')
    plt.ylim([bandwidth.min(), bandwidth.max()])
    plt.legend(loc=2,frameon = False)
    plt.show()

    phase_space = phonons.phase_space.flatten(order='C')

    # Plot phase space vs frequency
    print('\n')
    plt.figure()
    plt.plot(frequency[3:], phase_space[3:], '.', markersize=10)
    plt.ylabel ("Phase space", fontsize=25, fontweight='bold')
    plt.xlabel("$\\nu$ (THz)", fontsize=25, fontweight='bold')
    plt.show()

    # Direct access to properties
    # calculated during the simulation
    frequency = phonons.frequency.flatten(order='C')
    bandwidth = phonons.bandwidth.flatten(order='C')
    participation = phonons.participation_ratio.flatten(order='C')

    # Plot phonon bandwidth vs frequency
    print('\n')
    plt.plot(frequency[3:],bandwidth[3:],'.',markersize=10,label= 'broadening shape: ' + str(phonons.broadening_shape))
    plt.ylabel('$\Gamma$ (THz)', fontsize=25, fontweight='bold')
    plt.xlabel("$\\nu$ (THz)", fontsize=25, fontweight='bold')
    plt.ylim([bandwidth.min(), bandwidth.max()])
    plt.legend(loc=2,frameon = False)
    plt.show()

    phase_space = phonons.phase_space.flatten(order='C')

    # Plot phase space vs frequency
    print('\n')
    plt.figure()
    plt.plot(frequency[3:], phase_space[3:], '.', markersize=10)
    plt.ylabel ("Phase space", fontsize=25, fontweight='bold')
    plt.xlabel("$\\nu$ (THz)", fontsize=25, fontweight='bold')
    plt.show()

     # Plot phase space vs frequency
    print('\n')
    plt.figure()
    plt.plot(frequency[3:], participation[3:], '.', markersize=10)
    plt.ylabel ("PR", fontsize=25, fontweight='bold')
    plt.xlabel("$\\nu$ (THz)", fontsize=25, fontweight='bold')
    plt.show()

    from kaldo.conductivity import Conductivity

    # Calculate conductivity  with direct inversion approach (inverse)
    print('\n')
    inv_cond_matrix = (Conductivity(phonons=phonons, method='inverse').conductivity.sum(axis=0))
    print('Inverse conductivity (W/mK): %.3f'%(np.mean(np.diag(inv_cond_matrix))))
    print(inv_cond_matrix)

    # Calculate conductivity  with  relaxation time approximation (rta)
    print('\n')
    rta_cond_matrix = Conductivity(phonons=phonons, method='rta').conductivity.sum(axis=0)
    print('Rta conductivity (W/mK): %.3f'%(np.mean(np.diag(rta_cond_matrix))))
    print(rta_cond_matrix)
    # Calculate conductivity  with  self-consistent approach (sc)

    print('\n')
    sc_cond_matrix = Conductivity(phonons=phonons, method='sc',n_iterations=20).conductivity.sum(axis=0)
    print('Self-consistent conductivity (W/mK): %.3f'%(np.mean(np.diag(sc_cond_matrix))))
    print(sc_cond_matrix)

    # Load in scattering rate
    scattering_rate = np.load(
        data_folder +
        folder+'/'+k_label+'/300/quantum/bandwidth.npy', allow_pickle=True)

    # Derive lifetime, which is inverse of scattering rate
    life_time = scattering_rate ** (-1)

    # Denote lists to intake mean free path in each direction
    mean_free_path = []
    for i in range(3):
        mean_free_path.append(np.loadtxt(
        data_folder +
        folder+'/'+k_label+'/300/quantum/inverse/mean_free_path_' +
        str(i) + '.dat'))

    # Convert list to numpy array and compute norm for mean free path
    # Convert the unit from angstrom  to nanometer
    mean_free_path = np.array(mean_free_path).T
    mean_free_path_norm = np.linalg.norm(
        mean_free_path.reshape(-1, 3), axis=1) / 10.0

    # Plot observables in subplot
    figure(figsize=(12, 3))
    subplot(1,3, 1)
    set_fig_properties([gca()])
    scatter(frequency.flatten(order='C'),
            life_time, facecolor='w', edgecolor='r', s=10, marker='s')
    yscale('log')
    ylabel(r'$\tau \ (ps)$', fontsize=14)
    xlabel('Frequency (THz)', fontsize=14)

    subplot(1,3, 2)
    set_fig_properties([gca()])
    scatter(frequency.flatten(order='C'),
            scattering_rate, facecolor='w', edgecolor='r', s=10, marker='d')
    ylabel(r'$\Gamma \ (THz)$', fontsize=14)
    xlabel('Frequency (THz)', fontsize=14)

    subplot(1,3, 3)
    set_fig_properties([gca()])
    scatter(frequency.flatten(order='C'),
            mean_free_path_norm, facecolor='w', edgecolor='r', s=10, marker='8')
    ylabel(r'$\lambda \ (nm)$', fontsize=14)
    xlabel('Frequency (THz)', fontsize=14)
    yscale('log')
    ylim([1e-4, 1e2])
    subplots_adjust(wspace=0.33)
    show()

    # Denote zeros to intake kappa tensor
    kappa_tensor = np.zeros([mean_free_path.shape[0], 3, 3])
    for i in range(3):
        for j in range(3):
            kappa_tensor[:, i, j] = np.loadtxt(data_folder +
        folder+'/'+k_label+'/300/quantum/inverse/conductivity_' +
        str(i) + '_' + str(j) + '.dat')

    # Sum over the 0th dimension to recover 3-by-3 kappa matrix
    kappa_matrix = kappa_tensor.sum(axis=0)
    print("Bulk thermal conductivity: %.1f W m^-1 K^-1\n"
          %np.mean(np.diag(kappa_matrix)))
    print("kappa matrix: ")
    print(kappa_matrix)
    print('\n')

    # Compute kappa in per mode and cumulative representations
    kappa_per_mode = kappa_tensor.sum(axis=-1).sum(axis=1)
    freq_sorted, kappa_cum_wrt_freq = cumulative_cond_cal(
        frequency.flatten(order='C'), kappa_tensor)
    lambda_sorted, kappa_cum_wrt_lambda = cumulative_cond_cal(
        mean_free_path_norm, kappa_tensor)

    # Plot observables in subplot
    figure(figsize=(12, 3))
    subplot(1,3, 1)
    set_fig_properties([gca()])
    scatter(frequency.flatten(order='C'),
            100*kappa_per_mode/np.sum(kappa_per_mode), facecolor='w', edgecolor='r', s=10, marker='>')
    gca().axhline(y = 0, color='k', ls='--', lw=1)
    # ylabel(r'$\kappa_{per \ mode}\;\left(\frac{\rm{W}}{\rm{m}\cdot\rm{K}}\right)$',fontsize=14)
    ylabel(r'$\kappa_{per \ mode}\;\left(\%\right)$',fontsize=14)
    xlabel('Frequency (THz)', fontsize=14)
    # ylim([-0.5, 6])

    subplot(1,3, 2)
    set_fig_properties([gca()])
    # gca().axhline(y = 147, color='b', ls='--', lw=2,  label=r'$\kappa_{LDA + BTE} \approx 147 \ W \ m^{-1} \ K^{-1} $')
    # gca().axhline(y = 145, color='c', ls='--', lw=2,  label=r'$\kappa_{PBE + BTE} \approx 145 \ W \ m^{-1} \ K^{-1} $')
    plot(freq_sorted, kappa_cum_wrt_freq, 'r',
         label=r'$\kappa_{PACE + \kappa ALDo}$')
    # fill_between(freq_sorted, 130, 160, color='0.5',alpha=0.2,
    #              label=r'$\kappa_{exp}$')
    ylabel(r'$\kappa_{cumulative, \omega}\;\left(\frac{\rm{W}}{\rm{m}\cdot\rm{K}}\right)$',fontsize=14)
    xlabel('Frequency (THz)', fontsize=14)
    legend(loc=4, frameon=False, fontsize=6)

    subplot(1,3, 3)
    set_fig_properties([gca()])
    plot(lambda_sorted, kappa_cum_wrt_lambda, 'r')
    xlabel(r'$\lambda \ (nm)$', fontsize=14)
    ylabel(r'$\kappa_{cumulative, \lambda}\;\left(\frac{\rm{W}}{\rm{m}\cdot\rm{K}}\right)$',fontsize=14)
    xscale('log')
    xlim(0.1,1e2)
    subplots_adjust(wspace=0.33)
    show()

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))

#main
