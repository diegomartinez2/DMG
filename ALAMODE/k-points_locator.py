import spglib
from ase.io import read
from ase.dft.kpoints import sc_special_paths, sc_special_points
import numpy as np

def get_lattice_type_from_spglib(cell):
    """
    Usa spglib para obtener el grupo espacial.
    Puedes mapear el grupo espacial para inferir el tipo de red (fcc, bcc, hcp, etc.).
    """
    dataset = spglib.get_symmetry_dataset(cell)
    spg_symbol = dataset['international']
    print(f"Grupo espacial detectado por spglib: {spg_symbol}")

    # Mapa básico ejemplos comunes (amplia según casos que necesites)
    group_to_lattice = {
        'Fm-3m': 'fcc',
        'Fmmm': 'orthorhombic',
        'Pm-3m': 'sc',
        'Im-3m': 'bcc',
        'P6_3/mmc': 'hcp',
        # Agrega más grupos espaciales y su red correspondiente si quieres
    }

    # Buscar coincidencia parcial de grupo espacial para determinar lattice_type
    for key in group_to_lattice:
        if spg_symbol.startswith(key):
            return group_to_lattice[key]

    # Si no se reconoce, devolver None
    return None

def main():
    # Cambia 'POSCAR' por ruta a tu fichero POSCAR o data de LAMMPS
    atoms = read('POSCAR')
    #atoms = read('data.lammps')

    # Construir la tupla para spglib: (celda, posiciones fraccionales, números atómicos)
    cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.get_atomic_numbers())

    lattice_type = get_lattice_type_from_spglib(cell)
    if lattice_type is None:
        print("No se pudo detectar automáticamente el tipo de red para ASE.")
        return

    print(f"Tipo de lattice determinado: {lattice_type}")

    # Obtener rutas y puntos especiales para ese lattice_type
    if lattice_type not in sc_special_paths or lattice_type not in sc_special_points:
        print(f"No hay rutas estándar definidas en ase.dft.kpoints para '{lattice_type}'")
        return

    paths = sc_special_paths[lattice_type]
    points = sc_special_points[lattice_type]

    print("\nRutas de alta simetría disponibles:")
    for path in paths:
        print(' - ' + '-'.join(path))

    # Ejemplo: mostrar coordenadas de la primera ruta
    first_path = paths[0]
    print(f"\nCoordenadas fraccionales para la ruta: {'-'.join(first_path)}")
    for label in first_path:
        kpt = points[label]
        print(f"{label}: {kpt}")

    # Conversión simple a coordenadas absolutas reciprocales (en Å⁻¹)
    recip_cell = atoms.cell.reciprocal()
    print("\nCoordenadas cartesianas (recíprocas) de la primera ruta:")
    for label in first_path:
        kpt_frac = points[label]
        kpt_cart = np.dot(kpt_frac, recip_cell)
        print(f"{label}: {kpt_cart}")

#    print("Puntos especiales:")
#    for label, coord in points.items():
#        print(f"{label}: {coord}")

if __name__ == "__main__":
    main()
