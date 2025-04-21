#Example of LAMMPS calculation with ASE
#this code only outputs the potential energy.
#!/usr/local/bin/python
from ase import Atom, Atoms
from ase.build import bulk
from ase.calculators.lammpsrun import LAMMPS

parameters = {'pair_style': 'eam/alloy',
            'pair_coeff': ['* * NiAlH_jea.eam.alloy H Ni']}

files = ['NiAlH_jea.eam.alloy']

Ni = bulk('Ni', cubic=True)
H = Atom('H', position=Ni.cell.diagonal()/2)
NiH = Ni + H

lammps = LAMMPS(files=files, **parameters)

NiH.calc = lammps
print("Energy ", NiH.get_potential_energy())
