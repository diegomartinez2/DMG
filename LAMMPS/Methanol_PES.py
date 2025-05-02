#!/usr/local/bin/python
import os
import lammps
import numpy as np
import pandas as pd

"""
Here is an example on how to use LAMMPS in order to compute the energy of your system, considering a given methodology (molecular mechanics here) and a series of molecular structures.

The system consists in a methanol molecule in vacuum and the corresponding data file is methanol.data.

In the following input we simply read this data file and define the force field. As we will not do any simulation later, Verlet scheme and neighbor list options are not specified and default values are used. Of course, atom_style and force fields options must be consistent with the data file.
"""
methanol_inp = """# LAMMPS input file for Methanol

units real
atom_style full
boundary  f f f

bond_style harmonic
angle_style harmonic
dihedral_style opls
pair_style      lj/cut/coul/cut 14

read_data methanol.data

# LJ: OPLS
# energy = kcal/mol
# distance = angstrom
pair_coeff        1        1    0.0660     3.5000
pair_coeff        2        2    0.0300     2.5000
pair_coeff        3        3    0.1700     3.1200
pair_coeff        4        4    0.0000     0.0000
pair_modify shift no mix geometric
special_bonds amber

fix NVE all nve
"""


#One single calculation

#Again, we create the lammps object. Here we reduce output and redirect stdout to devnull.


lmp.close()
lmp = lammps.lammps(cmdargs=["-log", "none", "-screen", os.devnull,  "-nocite"])
lmp.commands_string(methanol_inp)

#Let’s make an energy calculation and output energy quantities.

lmp.command("run 0")

for energy in ["etotal", "evdwl", "ecoul", "ebond", "eangle", "edihed"]:
    print(f"{energy:8s} {lmp.get_thermo(energy):10.4f} kcal.mol-1")

"""It’s possible to extract the last (current) geometry of the calculation and reconstruct the molecule. Note that, in this particular case, with run 0 the geometry is the same as the input one.
Remember also, that in LAMMPS, by default, atoms are not sorted (if you perform several steps). It is thus safer to get atom ids.
Here we use the numpy property of the lammps object to get numpy arrays (more convenient) instead of pointers."""

coords = lmp.numpy.extract_atom("x")
print("coords\n", coords)

ids = lmp.numpy.extract_atom("id")
print("atom ids", ids)

i_types = lmp.numpy.extract_atom("type")
print("atom types", i_types)

elements = ["C", "H", "O", "H"]
dict_types = {i: el for i, el in enumerate(elements, 1)}

# TODO: sort by id if needed
print(f"\n{'el':2s} {'id':>4s}{'x':>12s} {'y':>12s} {'z':>12s}")
for iat in range(lmp.get_natoms()):
    specie = dict_types[i_types[iat]]
    xyz = coords[iat]
    line =  f"{specie:2s} {ids[iat]:4d}"
    line += " ".join([f"{x:12.6f}" for x in xyz])
    print(line)
lmp.close()

"""PES calculations
Hereafter, we compute the energy on a series of geometries read from an xyz file.
First we read the coordinates in the scan.allxyz file. It consists in 73 geometries of methanol along the dihedral angle H-O-C-H."""

scan_coords = list()
with open("scan.allxyz", "r") as f:
    chunks = f.read().split(">")
    for chunk in chunks:
        lines = chunk.strip().split("\n")
        coords = list()
        for line in lines[2:]:
            coords.append(line.split()[1:])
        scan_coords.append(coords)
scan_coords = np.array(scan_coords, dtype=np.float64)

n_geo, n_atoms, n_dim = scan_coords.shape
print(n_geo, n_atoms, n_dim)

#Now, again we initialize the LAMMPS object and load the force fields information as previously.

lmp.close()
lmp = lammps.lammps(cmdargs=["-log", "none", "-screen", os.devnull,  "-nocite"])
lmp.commands_string(methanol_inp)

"""We write a loop over the geometries read in the the scan.allxyz file. The atoms’ coordinates are accessible from a pointer by interacting with the current LAMMPS execution.

Hereafter, we get back the pointer in the p_coords variable. Then, in a loop, the coordinates are updated and the energy is computed as done previously for one single geometry.

Thanks to the default thermo variables, for each geometry we store the different energy components in a dictionary."""

# get the pointers on coordinates => used to update geometry
p_coords = lmp.extract_atom("x")

scan_data = list()
dihedral = 0
for step_coords in scan_coords:
    dihedral += 5.0

    # update coordinates
    # p_coords is not a numpy array but a double C pointer
    for iat in range(n_atoms):
        p_coords[iat][0], p_coords[iat][1], p_coords[iat][2] = step_coords[iat]

    # compute energy
    lmp.command("run 0")

    # get back results
    step_data = dict(dihedral=dihedral)
    for energy in ["etotal", "evdwl", "ecoul", "ebond", "eangle", "edihed"]:
        step_data[energy] = lmp.get_thermo(energy)

    scan_data.append(step_data)

lmp.close()

#We gather data in a pandas DataFrame in combination with QM energies associated to the geometries. The energies are available in the scan.qm.dat file.

df = pd.DataFrame(scan_data)
df.set_index("dihedral", inplace=True)
df.etotal -= df.etotal.min()
df.head()

qm_df = pd.read_csv(
    "scan.qm.dat",
    index_col=0,
    delim_whitespace=True,
    names=["dihedral", "QM"]
)
qm_df -= qm_df.QM.min()
qm_df *= 627.5095  # Hartree -> kcal.mol-1
qm_df.head()

"""Now depending on what you want to do, you can compare the contributions of each energy to the total energy and compare them to the QM energy.

Here, in the case of OPLS force field of a dihedral angle, considering the height of the energy barriers, differences are smaller than the expected accuracy."""

ax = df.plot(y=["edihed", "etotal"], marker="o", ls="--")
ax = qm_df.plot(ax=ax, marker="d", ls=":")
