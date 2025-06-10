#!/usr/bin/env python
from ase import Atoms
from ase.calculators.calculator import Calculator
import numpy as np

class CustomCalculator(Calculator):
    implemented_properties = ['energy', 'forces']  # List of properties your calculator supports

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def calculate(self, atoms, properties):
        self.atoms = atoms
        # Replace with your own calculation logic
        self.properties = {'energy': -100.0, 'forces': np.zeros((len(atoms), 3))}

    def get_potential_energy(self):
        return self.properties['energy']

    def get_forces(self):
        return self.properties['forces']

# Create an Atoms object
atoms = Atoms('He', positions=[[0, 0, 0]])

# Create an instance of your custom calculator
calculator = CustomCalculator()

# Attach the calculator to the Atoms object
atoms.calc = calculator

# Calculate properties
energy = atoms.get_potential_energy()
forces = atoms.get_forces()

print(f"Energy: {energy} eV")
print(f"Forces:\n{forces}")
