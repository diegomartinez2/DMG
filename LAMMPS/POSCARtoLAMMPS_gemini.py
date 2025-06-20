#!/usr/bin/env python
from ase.io import read, write

# Replace 'POSCAR' with the name of your VASP POSCAR file
poscar_file = 'POSCAR'
# Replace 'lammps.data' with the desired output LAMMPS data file name
lammps_data_file = 'lammps.data'
# Define atom types if they are not explicitly in your POSCAR (e.g., if you only have element symbols)
# If your POSCAR has element names, ASE usually handles it.
# Otherwise, you might need to map them:
# atom_type_mapping = {'Fe': 1, 'C': 2} # Example mapping: Fe is type 1, C is type 2

# Read the VASP POSCAR file
atoms = read(poscar_file, format='vasp')

# Assign atom types (important for LAMMPS).
# If your POSCAR already defines atom types (e.g., Fe C), ASE will usually
# assign them automatically as 1, 2, etc. based on their order.
# If you need specific integer types or want to ensure proper mapping:
# You might need to iterate through atoms and set atoms.info['atom_type']
# or ensure atoms.symbols correctly reflect the order.

# Write to LAMMPS data format
# 'atom_style' is crucial. Common choices are 'atomic', 'full', 'charge', 'molecular', etc.
# 'atomic' is basic (atom ID, type, x, y, z)
# 'full' includes molecule ID, atom ID, type, charge, x, y, z, etc.
# Choose based on your LAMMPS simulation needs.
write(lammps_data_file, atoms, format='lammps-data', atom_style='charge')

print(f"Successfully converted {poscar_file} to {lammps_data_file} with atom_style='charge'")
print("Remember to adjust atom_style and other parameters in the script if needed.")
print("You may also need to manually add atom masses to the LAMMPS data file if not automatically included or specify them in the 'write' function.")
