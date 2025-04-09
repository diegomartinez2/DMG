from ase import Atoms, Atom
from ase.calculators.lammpsrun import LAMMPS

a = [6.5, 6.5, 7.7]
d = 2.3608
NaCl = Atoms([Atom('Na', [0, 0, 0]),
                Atom('Cl', [0, 0, d])],
                cell=a, pbc=True)

calc = LAMMPS(tmp_dir="/tmp/lammpstmp")
NaCl.set_calculator(calc)

print(NaCl.get_stress())
