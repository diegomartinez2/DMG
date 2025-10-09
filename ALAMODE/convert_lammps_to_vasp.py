import sys
import os
import subprocess
from tempfile import NamedTemporaryFile
import numpy as np

# Note: This script assumes you have ASE (Atomic Simulation Environment) installed for the first method.
# Install with: pip install ase
# For the alternative method, install Open Babel: sudo apt install openbabel (on Ubuntu) or equivalent.
# The input LAMMPS file is read from stdin or a file path provided as argument.

def lammps_to_vasp_ase(lammps_file_path):
    """
    Method 1: Using ASE to convert LAMMPS data file to VASP POSCAR format.
    """
    try:
        from ase.io import read, write
        from ase import Atoms
    except ImportError:
        print("Error: ASE is not installed. Install with 'pip install ase'.")
        return None

    # Read the LAMMPS structure
    atoms = read(lammps_file_path, format='lammps-data')

    # Create a temporary file for VASP output
    with NamedTemporaryFile(mode='w', suffix='.vasp', delete=False) as tmp:
        vasp_file = tmp.name
        write(vasp_file, atoms, format='vasp', vasp5=True, sort=True, direct=True)

    # Read and print the content
    with open(vasp_file, 'r') as f:
        content = f.read()
    os.unlink(vasp_file)  # Clean up

    return content

def lammps_to_vasp_babel(lammps_file_path):
    """
    Method 2: Using Open Babel (via subprocess) to convert LAMMPS data to VASP POSCAR.
    Open Babel supports LAMMPS as 'lammpsdata' and VASP as 'poscar'.
    """
    try:
        # Create temporary files
        with NamedTemporaryFile(mode='w', suffix='.lammps', delete=False) as in_tmp:
            # Copy the LAMMPS content to temp file if needed, but here we assume input is file
            pass  # Input is already a file
            in_file = lammps_file_path

        out_file = NamedTemporaryFile(mode='w', suffix='.poscar', delete=False).name

        # Run Open Babel command
        cmd = [
            'obabel',
            '-ilammpsdata', in_file,
            '-oposcar', '-O', out_file
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Error running Open Babel: {result.stderr}")
            return None

        # Read and print the content
        with open(out_file, 'r') as f:
            content = f.read()
        os.unlink(out_file)  # Clean up

    except FileNotFoundError:
        print("Error: Open Babel is not installed or not in PATH. Install with 'sudo apt install openbabel'.")
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None

    return content

def main():
    if len(sys.argv) != 2:
        print("Usage: python convert_lammps_to_vasp.py <lammps_file_path>")
        sys.exit(1)

    lammps_file = sys.argv[1]

    if not os.path.exists(lammps_file):
        print(f"Error: File {lammps_file} not found.")
        sys.exit(1)

    print("=== Method 1: Using ASE ===")
    vasp_ase = lammps_to_vasp_ase(lammps_file)
    if vasp_ase:
        print(vasp_ase)

    print("\n=== Method 2: Using Open Babel ===")
    vasp_babel = lammps_to_vasp_babel(lammps_file)
    if vasp_babel:
        print(vasp_babel)

if __name__ == "__main__":
    main()
