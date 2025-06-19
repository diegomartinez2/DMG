#!/usr/bin/env python
# coding: utf-8

# Script to calculate the DOS and ... using kALDo, from a LAMMPS dynamical matrix


from ase.build import bulk, surface
from ase.build import cut
import ase.geometry
import numpy as np
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import BFGS
from ase.io import write, lammpsdata
from ase.visualize import view
from kaldo.forceconstants import ForceConstants
from ase.calculators.lammpslib import LAMMPSlib
from kaldo.phonons import Phonons
import kaldo.controllers.plotter as plotter
from pylab import *
import subprocess
from kaldo.conductivity import Conductivity
import os
#import lammps
import matplotlib.pyplot as plt

# In[2]:
#-------only for 2D materials--------------------
#Lz = 1      #this is the z of the calculation box in LAMMPS
#dist = 1    #this is the thickness of the 2D material
#print("Parameters for the conductivity calculation: Lz=",Lz," dist=",dist)
#------------------------------------------------
def cumulative_cond_cal(observables, kappa_tensor, prefactor=1/3):

    """Compute cumulative conductivity based on either frequency or mean-free path

       input:
       observables: (ndarray) either phonon frequency or mean-free path
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

# Denote plot default format
aw = 2
fs = 12
font = {'size': fs}
matplotlib.rc('font', **font)
matplotlib.rc('axes', linewidth=aw)

# Configure Matplotlib to use a LaTeX-like style without LaTeX
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'sans'
plt.rcParams['mathtext.fontset'] = 'cm'


# In[ ]:
#Supercell structure (this makes the supercell bigger for better calculations)

nx = 3
ny = 3
nz = 3
print("Supercell:",nx,"x",ny,"x",nz)

# In[19]:


#get_ipython().system('rm -r *fd* *xyz ge-* plots *lmp *log* graphene')

#atoms = lammpsdata.read_lammps_data('relax.dat',Z_of_type={1:6},sort_by_id=True,style='atomic',units='metal')
print("Read LAMMPS relaxed data")
atoms = lammpsdata.read_lammps_data('relax.dat',style='charge',units='metal')

supercell = np.array([nx,ny,nz])
print("Save replicated_atoms.xyz file")
write('replicated_atoms.xyz',format ='extxyz',images=atoms.copy().repeat(supercell))
print("set force constants from LAMMPS format data:")
print("from folder; FORMAT:lammps,Config_filename:replicated_atoms.xyz,ASE_format:xyz,2nd_order_FC:Dyn.form,3rd_order_FC:THIRD")
forceconstants = ForceConstants.from_folder(folder='./', supercell=supercell,format='lammps')


# In[21]:

print("Make a figure of the atoms positions")
from ase.visualize.plot import plot_atoms
fig, ax = plt.subplots()
plot_atoms(atoms.copy().repeat(supercell), ax, radii=0.3, rotation=('0x,0y,0z'))
fig.savefig("ase_atoms.png")
#savefig("ase_atoms.png")

# In[22]:

print("Set the kappa grid:")
kx=10
ky=10
kz=10
kpts = [kx, ky, kz]
temperature = 300   #WARNING: note that if this changes then you must change the sucesive directories in this code...
is_classic = False
k_label = str(kx) + '_' + str(ky) + '_' + str(kz)
print(k_label)

print("create phonon object")
# Create a phonon object
phonons = Phonons(forceconstants=forceconstants,
                kpts=kpts,
                is_classic=is_classic,
                temperature=temperature,
                folder='Supercell',
                storage='numpy')


# In[23]:

print("set kpath and plot dispersion")
kpath=atoms.cell.bandpath(path='GXVRGY',npoints=200,
special_points={'G': [0, 0, 0], 'X': [1/2, 0, 0], 'Y': [0, 1/2, 0],
'Z': [0, 0, 1/2], 'R': [1/2, 1/2, 1/2], 'T': [0, 1/2, 1/2], 'U': [1/2, 0, 1/2],
'V': [1/2, 0, 1/2]})
plotter.plot_dispersion(phonons,with_velocity=True,is_showing=True,manually_defined_path=kpath)
savefig("ase_dispersion.png")
print("plot DOS")
plotter.plot_dos(phonons,p_atoms=None,bandwidth=2,n_points=200,filename='dos',is_showing=False)
savefig("ase_DOS.png")
## Cargar constantes de fuerza de segundo orden desde Dyn.form
# La API de kALDo tiene una función para esto, por ejemplo:
#fc2 = ForceConstants.from_lammps(file_name='Dyn.form', format='shengbte', replicated_atoms_file='replicated_atoms.xyz')

# Si necesitas constantes de tercer orden (para conductividad térmica, por ejemplo)
# fc3 = ForceConstants.from_lammps_third_order(file_name='THIRD', replicated_atoms_file='replicated_atoms.xyz')
# Ejemplo conceptual con kALDo (la sintaxis exacta puede variar)
# Asegúrate de tener una clase de red definida para tu sistema
# Por ejemplo, si tienes una red de la clase Lattice de kALDo
# lattice = kaldo.Lattice(unit_cell_filename='your_structure.xyz')
#GAMMA=0,0,0  X=0.5,0,0  Y=0,0.5,0  Z=0,0,0.5  R=0.5,0.5,0.5  T=0,0.5,0.5  U=0.5,0,0.5  V=0.5,0.5,0
# Si ya tienes fc2 cargado
# dispersions = fc2.calculate_dispersions(q_points=your_q_point_path)
# dispersions.plot()
# In[24]:


print('\n')
print("calculate conductivity (by inverse method)")
inv_cond_matrix = (Conductivity(phonons=phonons, method='inverse').conductivity.sum(axis=0))
#print('Conductivity from inversion (W/m-K): %.3f' % ((Lz/(dist))*np.mean(np.diag(inv_cond_matrix[0:2,0:2]))))
#print((Lz/(dist))*inv_cond_matrix)
print('Conductivity from inversion (W/m-K): %.3f' % (np.mean(np.diag(inv_cond_matrix)))) #for 3D materials
print(inv_cond_matrix)
"""
This calculates the conductiviti from inversion, as it is a 2D system it only takes x and y.
The z is calculated from:
The problem with Supercell is that
BTE requires the system volume at the denominator

but the cell used for the calculations has a very large Lz because it
is in pbc and we don't want the graphene sheet interacting with its
image. So we have to multiply by a factor Lz/d where d is the z
thickness of a graphene monolayer which conventionally is chosen as
the distance between two layers (I might have deleted that value in
the notebook by mistake) and for the Airebo potential in the example I
sent you is 3.4A.
In the kALDo example there is no need to scale by this factor because
it calculates the thermal conductivity for bulk silicon.
Also, as you might notice, to average the thermal conductivity tensor
I'm taking the trace of the tensor excluding the zz value, while the
example from the kALDo page does not.

print('Conductivity from inversion (W/m-K): %.3f' %
((Lz/(dist))*np.mean(np.diag(inv_cond_matrix[0:2,0:2]))))

vs

print('Inverse conductivity (W/mK):
%.3f'%(np.mean(np.diag(inv_cond_matrix))))


"""
#print('Conductivity from inversion (W/m-K): %.3f' % np.mean(np.diag(inv_cond_matrix[0:2,0:2])) )
#print(inv_cond_matrix)


# In[25]:

print("Load in group velocity and frequency data")
data_folder = "Supercell/"+k_label
# Load in group velocity and frequency data
frequency =  np.load(
    data_folder + '/frequency.npy',
    allow_pickle=True)

cv =  np.load(
    data_folder + '/300/quantum/heat_capacity.npy',
    allow_pickle=True)

group_velocity = np.load(
    data_folder + '/velocity.npy')

phase_space = np.load(data_folder +
 '/300/quantum/_ps_and_gamma.npy',
                      allow_pickle=True)[:,0]

# Compute norm of group velocity
# Convert the unit from angstrom / picosecond to kilometer/ second
group_velcotiy_norm = np.linalg.norm(
    group_velocity.reshape(-1, 3), axis=1) / 10.0

# Plot observables in subplot
figure(figsize=(12, 3))
subplot(1,3, 1)
set_fig_properties([gca()])
scatter(frequency.flatten(order='C')[3:], 1e23*cv.flatten(order='C')[3:],
        facecolor='w', edgecolor='r', s=10, marker='8')
ylabel (r"$C_{v}$ ($10^{23}$ J/K)")
xlabel('Frequency (THz)', fontsize=14)
ylim(0.9*1e23*cv.flatten(order='C')[3:].min(), 1.05*1e23*cv.flatten(order='C')[3:].max())

subplot(1 ,3, 2)
set_fig_properties([gca()])
scatter(frequency.flatten(order='C'),
        group_velcotiy_norm, facecolor='w', edgecolor='r', s=10, marker='^')
xlabel('Frequency (THz)', fontsize=14)
ylabel(r'$|v| \ (\frac{km}{s})$', fontsize=14)

subplot(1 ,3, 3)
set_fig_properties([gca()])
scatter(frequency.flatten(order='C'),
        phase_space, facecolor='w', edgecolor='r', s=10, marker='o')
xlabel('Frequency (THz)', fontsize=14)
ylabel('Phase space', fontsize=14)
subplots_adjust(wspace=0.33)
show()
savefig("ase_multiplot.png")


# In[26]:


# Load in scattering rate
scattering_rate = np.load(
    data_folder +
    '/300/quantum/bandwidth.npy', allow_pickle=True)

# Derive lifetime, which is inverse of scattering rate
life_time = scattering_rate ** (-1)

# Denote lists to intake mean free path in each direction
mean_free_path = []
for i in range(3):
    mean_free_path.append(np.loadtxt(
    data_folder +
    '/300/quantum/inverse/mean_free_path_' +
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
xscale('log')
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
ylim([1e-1, 1e5])
subplots_adjust(wspace=0.33)
show()
savefig("ase_multiplot2.png")

# In[27]:


# Denote zeros to intake kappa tensor
kappa_tensor = np.zeros([mean_free_path.shape[0], 3, 3])
for i in range(3):
    for j in range(3):
        kappa_tensor[:, i, j] = np.loadtxt(data_folder +
    '/300/quantum/inverse/conductivity_' +
    str(i) + '_' + str(j) + '.dat')

# Sum over the 0th dimension to recover 3-by-3 kappa matrix
kappa_matrix = kappa_tensor.sum(axis=0)
print("Bulk thermal conductivity: %.1f W m^-1 K^-1\n"
      %((Lz/(dist))*np.mean(np.diag(kappa_matrix[0:2,0:2]))))
     # %((1/1)*np.mean(np.diag(kappa_matrix[0:2,0:2]))))
print("kappa matrix: ")
print((Lz/(dist))*kappa_matrix)
#print((1/1)*kappa_matrix)
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
        (Lz/(dist))*kappa_per_mode, facecolor='w', edgecolor='r', s=10, marker='>')
        #(1/1)*kappa_per_mode, facecolor='w', edgecolor='r', s=10, marker='>')
gca().axhline(y = 0, color='k', ls='--', lw=1)
ylabel(r'$\kappa_{per \ mode}\;\left(\frac{\rm{W}}{\rm{m}\cdot\rm{K}}\right)$',fontsize=14)
xlabel('Frequency (THz)', fontsize=14)
# ylim([-0.5, 6])

subplot(1,3, 2)
set_fig_properties([gca()])
plot(freq_sorted, (Lz/(dist))*kappa_cum_wrt_freq, 'r')
#plot(freq_sorted, (1/1)*kappa_cum_wrt_freq, 'r')
ylabel(r'$\kappa_{cumulative, \omega}\;\left(\frac{\rm{W}}{\rm{m}\cdot\rm{K}}\right)$',fontsize=14)
xlabel('Frequency (THz)', fontsize=14)
legend(loc=4, frameon=False, fontsize=6)

subplot(1,3, 3)
set_fig_properties([gca()])
plot(lambda_sorted, (Lz/(dist))*kappa_cum_wrt_lambda, 'r')
#plot(lambda_sorted, (1/1)*kappa_cum_wrt_lambda, 'r')
xlabel(r'$\lambda \ (nm)$', fontsize=14)
ylabel(r'$\kappa_{cumulative, \lambda}\;\left(\frac{\rm{W}}{\rm{m}\cdot\rm{K}}\right)$',fontsize=14)
xscale('log')
xlim(1e-1,)
subplots_adjust(wspace=0.33)
show()
savefig("ase_multiplot3.png")

# In[ ]:
