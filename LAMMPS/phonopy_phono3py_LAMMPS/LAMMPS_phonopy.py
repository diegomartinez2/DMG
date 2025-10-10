import numpy as np
import os
import subprocess
from ase.io import read, write
from ase import Atoms
from collections import defaultdict
from typing import Dict, List, Tuple

# =======================================================
# 1. CLASE CONTROLADORA DE LAMMPS Y CONVERSOR DE FUERZAS
# =======================================================

class LammpsForceCalculator:
    """
    Clase para calcular fuerzas usando LAMMPS con un potencial ML
    y convertir la salida al formato FORCE_SETS de Phonopy/Phono3py.
    """
    def __init__(self, lammps_exe: str, ml_potential_file: str, ml_pair_style: str = "mlip"):
        """
        Inicializa el controlador.
        :param lammps_exe: Ruta al ejecutable de LAMMPS (ej. 'lmp_ml_compiled').
        :param ml_potential_file: Nombre del archivo del potencial ML (ej. 'ml.model').
        :param ml_pair_style: Estilo de potencial ML usado en LAMMPS (ej. 'mlip', 'snap').
        """
        self.lammps_exe = lammps_exe
        self.ml_potential_file = ml_potential_file
        self.ml_pair_style = ml_pair_style
        self.atom_type_map: Dict[str, int] = {} # Mapeo de elemento a tipo de átomo de LAMMPS

    def _generate_lammps_input(self, data_file: str, dump_file: str):
        """Genera el script de entrada de LAMMPS para calcular fuerzas."""
        input_script = f"""
# LAMMPS input script for force calculation
units metal
atom_style atomic
boundary p p p

# Lee la estructura (convertida de SPOSCAR/POSCAR por ASE)
read_data {data_file}

# Definición del potencial ML. ASE necesita los tipos ya definidos en el data file.
# El comando pair_coeff es solo de ejemplo, ajustar a su potencial específico.
pair_style {self.ml_pair_style}
pair_coeff * * {self.ml_potential_file}

# Computa las fuerzas por átomo (fx, fy, fz)
compute f_all all property/atom fx fy fz

# DUMP: Escribe las fuerzas, IDs y tipos. Usar el formato 'atom' para tener id, type, x, y, z, fx, fy, fz
# Escribimos un solo frame (timestep 0) en el archivo dump
dump 1 all custom 1 {dump_file} id type x y z fx fy fz
dump_modify 1 sort id

# Ejecutar un solo paso para calcular las fuerzas y guardarlas en el dump file
run 0
        """
        return input_script.strip()

    def _convert_lammps_dump_to_forces(self, dump_file: str, num_atoms: int) -> np.ndarray:
        """
        Lee el archivo dump de LAMMPS y extrae las fuerzas.
        Asume que el dump file fue ordenado por ID de átomo.
        """
        try:
            # Leer todas las líneas
            with open(dump_file, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise RuntimeError(f"El archivo de salida de LAMMPS '{dump_file}' no se encontró.")

        # Buscar la línea 'ITEM: ATOMS id type x y z fx fy fz' (encabezado del dump)
        try:
            atoms_start_line = [i for i, line in enumerate(lines) if 'ITEM: ATOMS' in line][-1]
        except IndexError:
            raise RuntimeError("Formato de dump de LAMMPS no reconocido.")

        # Las fuerzas comienzan después del encabezado (línea 1: ID, línea 2: # de átomos, línea 3: caja, línea 4: encabezado de datos)
        # La cabecera del dump tiene 9 líneas
        force_lines = lines[atoms_start_line + 1:]

        # Parsear las fuerzas: las columnas 7, 8 y 9 son fx, fy, fz (índices 6, 7, 8 en base cero)
        forces = np.zeros((num_atoms, 3))

        # El comando 'dump_modify 1 sort id' en el input de LAMMPS asegura que
        # los IDs estén en orden ascendente (1, 2, 3, ...) lo que coincide con el
        # orden de átomos en la celda atómica/supercelda de Phonopy.
        for i, line in enumerate(force_lines):
            parts = line.split()
            if len(parts) < 9: continue # Ignorar líneas vacías o incompletas

            # parts[0] es el ID de átomo. Si está ordenado, i+1 == ID.
            # forces[i] = [fx, fy, fz]
            forces[i] = [float(parts[6]), float(parts[7]), float(parts[8])]

            if i >= num_atoms -1: break # Parar después de leer todos los átomos

        if forces.shape[0] != num_atoms:
            raise RuntimeError(f"Error en el número de átomos leídos. Esperado: {num_atoms}, Leído: {forces.shape[0]}.")

        return forces

    def calculate_forces_and_write_force_sets(self, displacement_directory: str,
                                              supercell_file: str = "SPOSCAR",
                                              displacement_vector: np.ndarray = None,
                                              is_fc3: bool = False):
        """
        Ejecuta el ciclo completo: Lee la estructura, calcula con LAMMPS y escribe FORCE_SETS.

        :param displacement_directory: Carpeta de trabajo (ej. 'disp-001' o 'disp_fc3-001').
        :param supercell_file: Nombre del archivo de estructura (ej. 'SPOSCAR' de Phonopy/Phono3py).
        :param displacement_vector: Vector de desplazamiento para el encabezado de FORCE_SETS.
        :param is_fc3: Booleano para indicar si escribir FORCE_SETS_3RD.
        """
        print(f"-> Procesando directorio: {displacement_directory}")

        # 1. Configuración de Archivos
        current_dir = os.getcwd()
        work_dir = os.path.join(current_dir, displacement_directory)

        # Nombres de archivos temporales
        lammps_data_file = "supercell.data"
        lammps_input_file = "in.force"
        lammps_dump_file = "forces.dump"
        force_set_file = "FORCE_SETS" if not is_fc3 else "FORCE_SETS_3RD"

        if not os.path.exists(work_dir):
            raise FileNotFoundError(f"Directorio de desplazamiento no encontrado: {displacement_directory}")

        os.chdir(work_dir)

        try:
            # 2. Leer la estructura (SPOSCAR) con ASE
            atoms = read(supercell_file, format='vasp')
            num_atoms = len(atoms)

            # Asignar tipos de átomo LAMMPS (1, 2, 3...) a cada elemento
            # Esto es crucial ya que LAMMPS usa tipos numéricos
            unique_elements = sorted(list(set(atoms.get_chemical_symbols())))
            self.atom_type_map = {sym: i + 1 for i, sym in enumerate(unique_elements)}

            # 3. Escribir el archivo de datos de LAMMPS (Data File)
            # ASE lo hace automáticamente, mapeando elementos a tipos de LAMMPS.
            write(lammps_data_file, atoms, format='lammps-data', specorder=unique_elements)

            # 4. Generar y Escribir el input de LAMMPS
            lammps_input_content = self._generate_lammps_input(lammps_data_file, lammps_dump_file)
            with open(lammps_input_file, 'w') as f:
                f.write(lammps_input_content)

            # 5. Ejecutar LAMMPS
            print(f"   Ejecutando {self.lammps_exe}...")
            # subprocess.run(..., check=True) asegura que el script falle si LAMMPS falla
            subprocess.run([self.lammps_exe, "-in", lammps_input_file], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # 6. Leer la salida de LAMMPS y extraer las fuerzas
            forces = self._convert_lammps_dump_to_forces(lammps_dump_file, num_atoms)

            # 7. Escribir el archivo FORCE_SETS
            self._write_force_sets(force_set_file, forces, displacement_vector, is_fc3)
            print(f"   ✅ Archivo {force_set_file} generado exitosamente.")

        except Exception as e:
            print(f"   ❌ Error en el directorio {displacement_directory}: {e}")
        finally:
            # Volver al directorio original
            os.chdir(current_dir)

    def _write_force_sets(self, filename: str, forces: np.ndarray,
                          displacement_vector: np.ndarray, is_fc3: bool):
        """Escribe las fuerzas calculadas en el formato FORCE_SETS/FORCE_SETS_3RD."""

        # En el flujo de trabajo de Phonopy/Phono3py, cada carpeta de desplazamiento
        # tiene un único desplazamiento, por lo que el archivo FORCE_SETS contiene
        # la información de ese único desplazamiento.

        # Lectura de la información del desplazamiento (atom index y vector)
        if displacement_vector is None:
             raise ValueError("El vector de desplazamiento (dx, dy, dz) debe ser proporcionado.")

        # Los archivos FORCE_SETS requieren el *índice* del átomo desplazado.
        # En el flujo de trabajo de Phonopy/Phono3py, esta información se lee
        # del archivo 'disp.yaml' o 'disp_fc3.yaml' en el directorio principal.
        # Aquí, asumimos que el vector de desplazamiento es [índice, dx, dy, dz].
        # NOTA: Simplificamos, el índice debe ser 1-based.
        atom_index = int(displacement_vector[0]) + 1
        disp_vec = displacement_vector[1:]

        num_atoms = forces.shape[0]

        with open(filename, 'w') as f:
            # Número de desplazamientos contenidos en este archivo (siempre 1 en este enfoque)
            f.write(f"1\n")

            # --- Encabezado del Desplazamiento 1 ---
            # 1. Índice del átomo desplazado (1-based)
            f.write(f"{atom_index}\n")
            # 2. Vector de desplazamiento (dx, dy, dz)
            f.write(f"{disp_vec[0]:.16f} {disp_vec[1]:.16f} {disp_vec[2]:.16f}\n")

            # 3. Fuerzas resultantes en la supercelda
            for i in range(num_atoms):
                f.write(f"{forces[i, 0]:.16f} {forces[i, 1]:.16f} {forces[i, 2]:.16f}\n")


# =======================================================
# 2. FUNCIÓN PARA COORDINAR LA EXTRACCIÓN DE FC
# =======================================================

def run_phonon_calculation(lammps_calculator: LammpsForceCalculator,
                           is_fc3: bool = False,
                           disp_yaml_file: str = "disp.yaml"):
    """
    Coordina la lectura de disp.yaml/disp_fc3.yaml y llama al calculador para
    cada directorio de desplazamiento.
    """

    # Este es el paso clave: leer los archivos de desplazamiento generados
    # por phonopy/phono3py para saber qué carpeta procesar y qué vector usar.
    try:
        from yaml import safe_load # Necesitas 'pyyaml'
    except ImportError:
        print("Error: El paquete 'pyyaml' es necesario para leer disp.yaml. Instálalo.")
        return

    print(f"\n--- Leyendo información de desplazamientos de {disp_yaml_file} ---")

    with open(disp_yaml_file, 'r') as f:
        disp_data = safe_load(f)

    # Identificar la clave de desplazamiento ('displacements' o 'displacements' en el formato FC3)
    disp_key = 'displacements'
    if is_fc3:
         disp_key = 'fc_undisplacements' # o 'displacements', depende de la versión y el formato exacto

    if disp_key not in disp_data:
        raise KeyError(f"No se encontró la clave '{disp_key}' en {disp_yaml_file}. Verificar formato.")

    for entry in disp_data[disp_key]:
        # 'entry' contiene el 'atom', 'direction', y 'displacement'
        # Ejemplo de formato (simplificado, esto varía):
        # atom: 0, direction: 0, displacement: 0.01045

        # En el archivo generado por Phonopy/Phono3py, la información del desplazamiento es:
        # - 'number': Índice del desplazamiento (0-based)
        # - 'displacement': El vector [dx, dy, dz]
        # - 'original_index': El índice del átomo desplazado (0-based)

        disp_index = entry['number']
        disp_dir = f"disp-{disp_index + 1:03d}" if not is_fc3 else f"disp_fc3-{disp_index + 1:03d}"

        # El vector de desplazamiento contiene la información necesaria para el encabezado
        # Vector: [original_index, dx, dy, dz] (donde dx, dy, dz es el desplazamiento REAL en esa estructura)
        # En la práctica, el código debe reconstruir el vector de desplazamiento
        # a partir de 'disp.yaml' o 'disp_fc3.yaml'.

        # Simplificamos: El archivo de desplazamiento de Phonopy/Phono3py ya contiene
        # todos los desplazamientos en la lista. Se debe procesar para obtener
        # el índice del átomo desplazado (original_index) y el vector de desplazamiento.

        # Reconstrucción del vector de desplazamiento (ejemplo ilustrativo):
        atom_index_0based = entry['atom']
        disp_vec_x = entry.get('displacement', [0.0, 0.0, 0.0])[0]
        disp_vec_y = entry.get('displacement', [0.0, 0.0, 0.0])[1]
        disp_vec_z = entry.get('displacement', [0.0, 0.0, 0.0])[2]

        # El archivo 'disp.yaml' es más complejo y define qué átomo se desplazó
        # y cuánto se movió. Se debe leer el vector final de desplazamiento.

        # Simplificación: Usamos el índice del átomo (0-based) y un vector de desplazamiento
        # Aquí pasamos un placeholder: [índice_0based, dx, dy, dz]
        # **ESTO DEBE SER AJUSTADO SEGÚN EL FORMATO REAL DE SUS ARCHIVOS YAML**
        disp_info = np.array([atom_index_0based, disp_vec_x, disp_vec_y, disp_vec_z])

        lammps_calculator.calculate_forces_and_write_force_sets(
            displacement_directory=disp_dir,
            supercell_file="SPOSCAR", # El archivo de estructura generada
            displacement_vector=disp_info,
            is_fc3=is_fc3
        )

# =======================================================
# 3. EJECUCIÓN PRINCIPAL
# =======================================================

if __name__ == "__main__":
    # --- CONFIGURACIÓN DEL USUARIO ---
    LAMMPS_EXE = "lmp_ml_compiled"  # Asegúrate de que este ejecutable sea accesible
    ML_POTENTIAL = "ml.model"       # Archivo de potencial ML
    ML_STYLE = "mlip"               # Estilo de potencial (mlip, snap, pace, etc.)

    # -----------------------------------

    calculator = LammpsForceCalculator(
        lammps_exe=LAMMPS_EXE,
        ml_potential_file=ML_POTENTIAL,
        ml_pair_style=ML_STYLE
    )

    # Nota: Este script debe ejecutarse *después* de que Phonopy/Phono3py hayan generado
    # las carpetas de desplazamiento (disp-* / disp_fc3-*) y los archivos YAML.

    # 1. Cálculo de FC2 (Armónico)
    try:
        run_phonon_calculation(
            lammps_calculator=calculator,
            is_fc3=False,
            disp_yaml_file="disp.yaml"
        )
    except Exception as e:
        print(f"Fallo en el cálculo de FC2: {e}")

    # 2. Cálculo de FC3 (Anarmónico)
    try:
        run_phonon_calculation(
            lammps_calculator=calculator,
            is_fc3=True,
            disp_yaml_file="disp_fc3.yaml" # Nombre típico del archivo FC3
        )
    except Exception as e:
        print(f"Fallo en el cálculo de FC3: {e}")
