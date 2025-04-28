#Example of LAMMPS calculation with ASE
#this code only outputs the potential energy.
"""the code now tries to output the total energy and can output:
etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press
"""
#!/usr/local/bin/python
from ase import Atom, Atoms
from ase.build import bulk
from ase.calculators.lammpsrun import LAMMPS

parameters = {'pair_style': 'eam/alloy',
            'pair_coeff': ['* * NiAlH_jea.eam.alloy H Ni']}

files = ['NiAlH_jea.eam.alloy'] #included in lammps/potentials as examples
#The prefix of each file indicates the element(s)
#The suffix of each file indicates the pair style it is used

Ni = bulk('Ni', cubic=True)
H = Atom('H', position=Ni.cell.diagonal()/2)
NiH = Ni + H

lammps = LAMMPS(files=files, **parameters)

NiH.calc = lammps #what kind of calculation?, NVE? (calculation with default settings?)
print("Energy ", NiH.get_potential_energy())
input_lammps = """
thermo_style multi #Style multi prints a multiple-line listing of thermodynamic info that is the equivalent of “thermo_style custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press”. The listing contains numeric values and a string ID for each quantity.
thermo 200
"""
#lammps.commands_string(input_lammps) #testing adding termo
print("Total energy", NiH.get_thermo("etotal")) #does this work without setting the thermo_style?
