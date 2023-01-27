from ase.io import read, write

# outcar contanining one (single-point) or multiple (trajectory) DFT frames
in_filename = './OUTCAR'
out_filename = 'nequip-data.extxyz'

# read all frames into a list of ase.Atoms objects
all_atoms = read(in_filename, format='vasp-out', index=':')

for curr_atoms in enumerate(all_atoms): 
  write('out_filename', curr_atoms, append=True, format='extxyz')
