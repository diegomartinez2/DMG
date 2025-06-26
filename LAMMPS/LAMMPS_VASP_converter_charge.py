#!/usr/bin/env python
"""
This code uses ASE to read and write data files of atomic positions.
It can convert VASP (POSCAR) to LAMMPS data files (atomic or charge style).
It can convert LAMMPS data files (atomic or charge style) back to VASP (POSCAR).

Requires the 'ase' library: pip install ase

Note on LAMMPS 'charge' style:
- When converting from POSCAR to LAMMPS 'charge' style, POSCAR files typically
  do NOT contain charge information. By default, this script will assign a charge
  of 0.0 to all atoms. You can provide a list of custom charges if needed.
- When converting from LAMMPS 'charge' style to POSCAR, the charge information
  will be read by ASE but will NOT be written to the POSCAR file, as POSCAR
  does not natively support storing charges.
"""
import numpy as np
from ase.io import read, write
from ase import Atoms
from collections import OrderedDict

# Dictionary of atomic masses (in amu) for common elements
# ASE usually handles masses automatically from chemical symbols,
# but this dictionary can be used as a fallback or for verification.
ATOMIC_MASSES = {
    'H': 1.00794, 'He': 4.00260, 'Li': 6.941, 'Be': 9.01218, 'B': 10.811,
    'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'F': 18.9984, 'Ne': 20.1797,
    'Na': 22.9897, 'Mg': 24.305, 'Al': 26.9815, 'Si': 28.0855, 'P': 30.9738,
    'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
    'Sc': 44.9559, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938,
    'Fe': 55.845, 'Co': 58.9332, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38,
    'Ge': 72.64, 'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.80,
    'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.9059, 'Zr': 91.224, 'Nb': 92.9064,
    'Mo': 95.96, 'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
    'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.710, 'Sb': 121.760,
    'I': 126.9045, 'Xe': 131.293, 'Cs': 132.9055, 'Ba': 137.327, 'La': 138.9055,
    'Ce': 140.116, 'Pr': 140.9077, 'Nd': 144.242, 'Pm': 145.0, 'Sm': 150.36,
    'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.9253, 'Dy': 162.500, 'Ho': 164.9303,
    'Er': 167.259, 'Tm': 168.9342, 'Yb': 173.04, 'Lu': 174.967, 'Hf': 178.49,
    'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
    'Pt': 195.078, 'Au': 196.9665, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2,
    'Bi': 208.9804, 'Po': 209.0, 'At': 210.0, 'Rn': 222.0
}


def poscar_to_lammps(input_file, output_file, atom_style='atomic', charges=None):
    """
    Convert VASP POSCAR to LAMMPS data file using ASE.
    Supports 'atomic' and 'charge' atom_styles.

    Args:
        input_file (str): Path to the input VASP POSCAR file.
        output_file (str): Path to the output LAMMPS data file.
        atom_style (str): LAMMPS atom_style for the output file ('atomic' or 'charge').
                          Defaults to 'atomic'.
        charges (list, optional): List of charges for each atom. Required if atom_style
                                  is 'charge' and not all charges are 0.0.
                                  If None and atom_style is 'charge', all charges default to 0.0.
    """
    if atom_style not in ['atomic', 'charge']:
        raise ValueError("Unsupported atom_style. Choose 'atomic' or 'charge'.")

    # Read POSCAR file using ASE
    atoms = read(input_file, format='vasp')

    # Assign charges if 'charge' style is requested
    if atom_style == 'charge':
        if charges is None:
            # Default to zero charges if not provided
            atoms.set_initial_charges(np.zeros(len(atoms)))
            print(f"Warning: No charges provided for '{input_file}'. Assigning 0.0 charge to all atoms for LAMMPS 'charge' style.")
        elif len(charges) != len(atoms):
            raise ValueError(f"Number of charges ({len(charges)}) must match number of atoms ({len(atoms)}) in {input_file}.")
        else:
            atoms.set_initial_charges(charges)

    # Write to LAMMPS data file using ASE
    # ASE's lammps-data writer automatically handles Masses section from atom symbols
    # and also handles cell dimensions including triclinic cells.
    write(output_file, atoms, format='lammps-data', atom_style=atom_style)

    print(f"Converted {input_file} to {output_file} (LAMMPS '{atom_style}' style)")


def lammps_to_poscar(input_file, output_file, atom_symbols):
    """
    Convert LAMMPS data file (atomic or charge style) to VASP POSCAR using ASE.

    Args:
        input_file (str): Path to the input LAMMPS data file.
        output_file (str): Path to the output VASP POSCAR file.
        atom_symbols (list or dict): Mapping from LAMMPS atom types (1-based integer)
                                     to chemical symbols (e.g., ['Si', 'O'] or {1: 'Si', 2: 'O'}).
                                     This is crucial as LAMMPS data files only store integer types.
    """
    # Read LAMMPS data file using ASE.
    # ASE will automatically detect the atom_style ('atomic', 'charge', etc.)
    # It will store atom types in atoms.numbers and charges in atoms.arrays['charges'] if present.
    atoms = read(input_file, format='lammps-data')

    # Convert atom types (integers) back to chemical symbols
    # We need to create the 'symbols' attribute based on the provided mapping.
    if isinstance(atom_symbols, list):
        # If list, assume it's ordered by type_id-1 (e.g., atom_symbols[0] for type 1)
        if max(atoms.numbers) > len(atom_symbols) or min(atoms.numbers) < 1:
            raise ValueError("LAMMPS atom types are out of range for provided atom_symbols list.")
        mapped_symbols = [atom_symbols[t - 1] for t in atoms.numbers]
    elif isinstance(atom_symbols, dict):
        # If dict, map directly by type_id
        mapped_symbols = [atom_symbols[t] for t in atoms.numbers]
        if not all(t in atom_symbols for t in atoms.numbers):
            raise ValueError("Not all LAMMPS atom types found in provided atom_symbols dictionary.")
    else:
        raise TypeError("atom_symbols must be a list (ordered by type) or a dictionary (mapping type_id to symbol).")

    atoms.set_chemical_symbols(mapped_symbols)

    # Write POSCAR file using ASE.
    # Note: VASP POSCAR does not store charges directly.
    write(output_file, atoms, format='vasp', vasp5=True, direct=True)

    print(f"Converted {input_file} to {output_file} (VASP POSCAR).")
    if 'charges' in atoms.arrays:
        print("Note: Charge information was present in the LAMMPS data file but is NOT written to POSCAR.")


# --- Example Usage ---
if __name__ == "__main__":
    print("--- VASP POSCAR to LAMMPS Data Conversion Examples ---")

    # Example 1: Convert POSCAR to LAMMPS 'atomic' style
    # Create a dummy POSCAR file for demonstration if it doesn't exist
    poscar_content_si = """Silicon
1.0
5.4300000000000000   0.0000000000000000   0.0000000000000000
0.0000000000000000   5.4300000000000000   0.0000000000000000
0.0000000000000000   0.0000000000000000   5.4300000000000000
Si
8
Direct
0.0000000000000000  0.0000000000000000  0.0000000000000000
0.5000000000000000  0.5000000000000000  0.0000000000000000
0.5000000000000000  0.0000000000000000  0.5000000000000000
0.0000000000000000  0.5000000000000000  0.5000000000000000
0.2500000000000000  0.2500000000000000  0.2500000000000000
0.7500000000000000  0.7500000000000000  0.2500000000000000
0.7500000000000000  0.2500000000000000  0.7500000000000000
0.2500000000000000  0.7500000000000000  0.7500000000000000
"""
    with open('POSCAR_Si', 'w') as f:
        f.write(poscar_content_si)

    poscar_to_lammps('POSCAR_Si', 'data.lammps_atomic_Si', atom_style='atomic')

    # Example 2: Convert POSCAR (with two elements) to LAMMPS 'charge' style
    # Create a dummy POSCAR file for demonstration
    poscar_content_sio2 = """SiO2
1.0
4.9130000000000000   0.0000000000000000   0.0000000000000000
0.0000000000000000   8.4890000000000000   0.0000000000000000
0.0000000000000000   0.0000000000000000   5.4020000000000000
Si O
1 2
Direct
0.0000000000000000  0.0000000000000000  0.0000000000000000
0.5000000000000000  0.5000000000000000  0.5000000000000000
0.2500000000000000  0.2500000000000000  0.0000000000000000
0.7500000000000000  0.7500000000000000  0.0000000000000000
0.2500000000000000  0.7500000000000000  0.5000000000000000
0.7500000000000000  0.2500000000000000  0.5000000000000000
"""
    with open('POSCAR_SiO2', 'w') as f:
        f.write(poscar_content_sio2)

    # Let's assume charges for SiO2 (Si: +2.0, O: -1.0)
    # Number of Si atoms: 1, Number of O atoms: 2
    # So charges = [Si_charge, O_charge, O_charge]
    # For the provided POSCAR, it's 1 Si and 2 O, so total 3 atoms.
    # We need to know the order of atoms in the POSCAR.
    # ASE's read('vasp') will typically order by element, then by appearance.
    # The dummy POSCAR has Si first, then O. So the charges list should match:
    # 1 Si, 2 O => charges for 1st Si, 1st O, 2nd O...
    # Based on the dummy POSCAR, it's 1 Si, 2 O:
    # Si atom 1
    # O atom 2
    # O atom 3
    # A more robust way to assign charges by type would be to iterate `atoms.symbols` after reading the POSCAR
    # Example for POSCAR_SiO2 with 1 Si and 2 O:
    # Si @ (0,0,0) -> Type 1
    # O @ (0.5,0.5,0.5) -> Type 2
    # O @ (0.25,0.25,0.0) -> Type 2
    # The current ASE read for VASP will assign 1 Si and 2 O atoms
    # charges_for_sio2 = [2.0, -1.0, -1.0] # This assumes order Si, O, O
    # Let's generate a slightly more complex POSCAR for charging example to be clear.
    poscar_content_complex = """Mixed Oxide Example
1.0
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
Al O
2 3
Direct
0.0 0.0 0.0
0.5 0.5 0.5
0.25 0.25 0.0
0.75 0.25 0.0
0.25 0.75 0.0
"""
    with open('POSCAR_Mixed', 'w') as f:
        f.write(poscar_content_complex)

    # For POSCAR_Mixed: 2 Al atoms, 3 O atoms (total 5 atoms)
    # Default charges (e.g., Al +3.0, O -2.0)
    # The order of symbols in POSCAR_Mixed is Al, O.
    # So atoms will be Al, Al, O, O, O in that order from ASE read.
    charges_mixed_example = [3.0, 3.0, -2.0, -2.0, -2.0]
    poscar_to_lammps('POSCAR_Mixed', 'data.lammps_charge_Mixed', atom_style='charge', charges=charges_mixed_example)

    # Example 3: Convert POSCAR to LAMMPS 'charge' style with default zero charges
    poscar_to_lammps('POSCAR_Si', 'data.lammps_charge_Si_default_zero', atom_style='charge')


    print("\n--- LAMMPS Data to VASP POSCAR Conversion Examples ---")

    # Example 1: Convert LAMMPS 'atomic' data back to POSCAR
    # Need to provide the mapping from LAMMPS atom type (integer) to chemical symbol.
    # For data.lammps_atomic_Si, type 1 corresponds to 'Si'.
    lammps_atomic_symbols = {1: 'Si'} # Using dictionary for clarity
    lammps_to_poscar('data.lammps_atomic_Si', 'POSCAR_Si_converted_atomic', lammps_atomic_symbols)

    # Example 2: Convert LAMMPS 'charge' data back to POSCAR
    # For data.lammps_charge_Mixed, type 1 is 'Al', type 2 is 'O'.
    lammps_charge_symbols_mixed = {1: 'Al', 2: 'O'}
    lammps_to_poscar('data.lammps_charge_Mixed', 'POSCAR_Mixed_converted_charge', lammps_charge_symbols_mixed)

    # Example 3: Convert LAMMPS 'charge' data (with default zero charges) back to POSCAR
    lammps_charge_symbols_si_default = {1: 'Si'}
    lammps_to_poscar('data.lammps_charge_Si_default_zero', 'POSCAR_Si_converted_charge_default_zero', lammps_charge_symbols_si_default)

    print("\nConversion examples completed. Check the generated 'data.lammps_*' and 'POSCAR_*_converted' files.")
