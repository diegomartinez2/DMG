from ase.io import read, write
from ase import Atoms
import numpy as np

def poscar_to_lammps(poscar_file, output_file, atom_types=None):
    """
    Convert a POSCAR file to a LAMMPS data file using ASE.
    
    Parameters:
    poscar_file (str): Path to the POSCAR file.
    output_file (str): Path to the output LAMMPS data file.
    atom_types (dict): Dictionary mapping chemical symbols to integer atom types.
                       If None, types are assigned automatically in order of appearance.
    """
    # Read the POSCAR file
    atoms = read(poscar_file, format='vasp')

    # Get chemical symbols
    symbols = atoms.get_chemical_symbols()
    unique_symbols = list(dict.fromkeys(symbols))  # Preserve order of appearance

    # Assign atom types if not provided
    if atom_types is None:
        atom_types = {symbol: idx + 1 for idx, symbol in enumerate(unique_symbols)}
    
    # Set atom types as tags
    tags = [atom_types[symbol] for symbol in symbols]
    atoms.set_tags(tags)

    # Write to LAMMPS data file
    write(output_file, atoms, format='lammps-data')

    print(f"LAMMPS data file written to {output_file}")
    print("Atom types assigned:", atom_types)

# Example usage
if __name__ == "__main__":
    # Specify the input POSCAR file and output LAMMPS file
    poscar_file = 'POSCAR'
    output_file = 'structure.lammps'

    # Define atom types (optional, maps chemical symbols to LAMMPS atom types)
    atom_types = {'S': 1, 'C': 2, 'Cu': 3}

    # Convert POSCAR to LAMMPS
    poscar_to_lammps(poscar_file, output_file, atom_types=atom_types)