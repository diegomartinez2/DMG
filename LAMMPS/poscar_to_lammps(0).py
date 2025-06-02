#!/usr/local/bin/python
"""
This code tries to transform from vasp POSCAR to LAMMPS datafile.
It laks the masses and charges.
"""
from ase.io import read, write
from ase import Atoms
import numpy as np

def poscar_to_lammps(poscar_file, output_file, atom_types=None):
    """
    Convert a POSCAR file to a LAMMPS data file using ASE.

    Parameters:
    poscar_file (str): Path to the POSCAR file or string containing POSCAR content.
    output_file (str): Path to the output LAMMPS data file.
    atom_types (dict): Dictionary mapping chemical symbols to integer atom types.
                       If None, types are assigned automatically in order of appearance.
    """
    # Read the POSCAR file or string
    try:
        # If poscar_file is a file path, read it directly
        atoms = read(poscar_file, format='vasp')
    except FileNotFoundError:
        # If poscar_file is a string, write it to a temporary file
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.POSCAR', delete=False) as tmp:
            tmp.write(poscar_file)
            tmp_path = tmp.name
        atoms = read(tmp_path, format='vasp')
        import os
        os.remove(tmp_path)  # Clean up temporary file

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
    # POSCAR content as a string (from your input)
    poscar_content = """Cu3BHT Stacking AA
   1.0000000000000000
     8.6583927714873727   -0.0012334306861130   -0.0977516149958681
    -4.3297941313395762    7.5169175210519228    0.3625204288635386
     0.0199674387339069    1.0622118669401301    3.2308891839373146
   S    C    Cu
     6     6     3
Direct
  0.4182766067522778  0.2123490912898226  0.9693989501296163
  0.2139986826536781  0.4277357740577939  0.9199623281519470
  0.5817234222477208  0.7876508937101726  0.0306010468703797
  0.7860013023463170  0.5722642549422051  0.0800376688480491
  0.2037152632089848  0.7833079775212135  0.0632300419137548
  0.7962847517910130  0.2166920374787844  0.9367699360862506
  0.1882998155346910  0.0949355494124820  0.9941841082553077
  0.0960067921368322  0.1919451023892282  0.9761969349882244
  0.8117001984653102  0.9050644655875160  0.0058158717446979
  0.9039932228631656  0.8080549116107731  0.0238030450117885
  0.0925260714876669  0.9032551607408786  0.0190606584572488
  0.9074739215123359  0.0967448322591241  0.9809393275427570
  0.5000000000000000  0.0000000000000000 -0.0000000000000000
 -0.0000000000000000  0.5000000000000000 -0.0000000000000000
  0.5000000000000000  0.5000000000000000 -0.0000000000000000"""

    # Define atom types (optional, maps chemical symbols to LAMMPS atom types)
    atom_types = {'S': 1, 'C': 2, 'Cu': 3}

    # Convert POSCAR to LAMMPS
    poscar_to_lammps(poscar_content, 'structure.lammps', atom_types=atom_types)
