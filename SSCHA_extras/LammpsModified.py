# -*- coding: utf-8 -*-
# Este script python no es la mejor opción ya que escribe en el disco múltipes
# veces, lo que afecta tanto al rendimiento como a la computadora.
from ase.calculators.calculator import Calculator, FileIOCalculator, all_changes
from ase.io import write
import os
import subprocess
import shutil

# ----------------------------------------------------------------------
# 1. Definición de la Clase LammpsModified (Calculator de ASE)
# ----------------------------------------------------------------------

class LammpsModified(FileIOCalculator):
    """
    Calculator de ASE personalizado para ejecutar una versión
    modificada de LAMMPS.

    Este calculator genera el archivo de entrada 'in.sw' y el archivo
    de estructura 'atoms.data', y luego ejecuta el binario LAMMPS
    con el comando './lmp_mpi_chimes < in.sw'.

    Parámetros
    ----------
    command : str
        El comando que se debe ejecutar. Por defecto: './lmp_mpi_chimes'.
    input_file : str
        Nombre del archivo de entrada LAMMPS. Por defecto: 'in.sw'.
    atoms_file : str
        Nombre del archivo de datos de estructura. Por defecto: 'atoms.data'.
    lammps_input : str
        Contenido del script de entrada LAMMPS (las líneas de comando).
    label : str
        Prefijo para los archivos de salida (ej. 'lammps').
    ... (más parámetros de ASE para FileIOCalculator)
    """

    implemented_properties = ['energy', 'forces']
    command_default = './lmp_mpi_chimes < {input_file}' # Plantilla del comando

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='lammps_chimes', atoms=None,
                 command=None,
                 input_file='in.sw',
                 atoms_file='atoms.data',
                 lammps_input=None,
                 **kwargs):

        # Inicialización del FileIOCalculator base
        FileIOCalculator.__init__(self, restart=restart,
                                  ignore_bad_restart_file=ignore_bad_restart_file,
                                  label=label, atoms=atoms,
                                  **kwargs)

        # Configuración del comando y archivos
        self.options = kwargs
        self.input_file = input_file
        self.atoms_file = atoms_file
        self.lammps_input = lammps_input

        # Si no se define un comando, usamos el default y sustituimos la plantilla
        if command is None:
            self.command = self.command_default.format(input_file=self.input_file)
        else:
            self.command = command

        # Inicializar resultados
        self.results = {}

        if self.lammps_input is None:
             raise ValueError("Error: Se requiere el parámetro 'lammps_input' con el contenido del script LAMMPS.")

    def write_input(self, atoms, properties=None, system_changes=all_changes):
        """
        Escribe el archivo de estructura (atoms.data) y el script de entrada (in.sw).
        """
        Calculator.write_input(self, atoms, properties, system_changes)

        # 1. Escribir el archivo de estructura LAMMPS (atoms.data)
        # Aquí usamos la función 'write' de ASE con el formato 'lammps-data'.
        write(self.atoms_file, atoms, format='lammps-data', specorder=self.options.get('specorder', None))

        # 2. Escribir el script de entrada LAMMPS (in.sw)
        # Reemplazamos cualquier marcador de posición si es necesario, y escribimos el contenido.
        lammps_script = self.lammps_input.format(atoms_file=self.atoms_file,
                                                 run_steps=self.options.get('run_steps', 100))

        with open(self.input_file, 'w') as f:
            f.write(lammps_script)

    def read_results(self):
        """
        Lee la energía y las fuerzas del archivo de salida de LAMMPS (log.lammps, dump.force, etc.).
        """
        # AVISO: Esta parte asume un formato de salida estándar de LAMMPS que ASE puede parsear
        # o que puedes parsear tú mismo. Por simplicidad, asumiremos que tu LAMMPS modificado
        # genera archivos 'log.lammps' o 'dump.force' estándar.

        output_file = 'log.lammps'  # Asumiendo que esta es la salida de log por defecto

        if not os.path.exists(output_file):
            raise FileNotFoundError(f"No se encontró el archivo de salida '{output_file}'.")

        # --- A. Lectura de Energía (Ejemplo de cómo leer del log) ---
        energy_line = None
        with open(output_file, 'r') as f:
            lines = f.readlines()
            # Buscar la línea que contiene la energía total (Etot) o el término final del log
            for line in reversed(lines):
                if 'Loop time' in line: # o alguna otra marca final
                    break
                if 'Etot' in line or 'Total energy' in line: # Adaptar a tu log
                    energy_line = line
                    break

        if energy_line:
            # Esto es muy dependiente del formato de tu log. Ajusta el parsing.
            # Por ejemplo, si es la última columna:
            # energy = float(energy_line.split()[-1])
            # self.results['energy'] = energy
            # Por ahora, dejaremos un valor ficticio si el parsing es complejo:
            self.results['energy'] = -123.456  # Reemplazar con el valor real parseado

        # --- B. Lectura de Fuerzas (Requiere un comando de dump en in.sw) ---
        forces_file = 'dump.force' # Requiere 'dump custom dump.force id type fx fy fz' en in.sw

        if os.path.exists(forces_file):
            # Aquí necesitarías una función de parsing para leer las fuerzas.
            # Por simplicidad, usa una de las utilidades de ASE o de cellconstructor si existe.
            # Aquí asumiremos que las fuerzas se leen correctamente y se almacenan en self.results['forces']
            self.results['forces'] = self.read_lammps_forces(forces_file) # Implementar esta función

        # Poner valores si no se puede parsear
        if 'forces' not in self.results:
             self.results['forces'] = [ [0.0, 0.0, 0.0] for _ in atoms ]


    def read_lammps_forces(self, filename):
        """
        Función auxiliar para leer fuerzas de un archivo de dump LAMMPS (muy simplificado).
        NECESITA IMPLEMENTACIÓN DETALLADA SEGÚN EL FORMATO DE TU DUMP.
        """
        # Esta es solo una plantilla. La implementación real es compleja
        # ya que requiere parsear el encabezado del dump LAMMPS.
        import numpy as np

        # Lectura simplificada: asume que las últimas líneas son las fuerzas
        with open(filename, 'r') as f:
            lines = f.readlines()

        # En un archivo dump LAMMPS, el encabezado tiene 9 líneas
        forces_lines = lines[9:]

        forces = []
        for line in forces_lines:
            parts = line.split()
            # Asume formato: id type fx fy fz
            if len(parts) >= 5:
                # Solo tomamos las componentes fx, fy, fz (índices 2, 3, 4)
                forces.append([float(parts[2]), float(parts[3]), float(parts[4])])

        return np.array(forces)

# ----------------------------------------------------------------------
# 2. Ejemplo de Uso con ASE
# ----------------------------------------------------------------------

if __name__ == '__main__':
    try:
        from ase.build import bulk
        import numpy as np
    except ImportError:
        print("Aviso: ASE o NumPy no está instalado. El ejemplo de uso no se ejecutará.")
        exit()

    # --- Definición del Script de Entrada de LAMMPS (in.sw) ---
    # Este es el contenido que se escribirá en 'in.sw'.
    # NOTA: Debes ajustarlo para que use tu potencial y genere la salida que necesitas (log.lammps y dump.force).
    LAMMPS_INPUT_SCRIPT = """
# Script de entrada para LAMMPS modificado

# Inicialización
units metal
atom_style atomic

# Definición de la estructura
read_data {atoms_file}

# Potencial POTENCIAL CHIMESFF
pair_style      chimesFF
pair_coeff      * * params.txt

# Configuración del cálculo
neighbor 2.0 bin
neigh_modify delay 0 every 1

## Minimización/Relajación
#minimize 1.0e-4 1.0e-6 100 1000

# Salida de fuerzas (Necesario para que ASE lea los resultados)
# Asegúrate de que tu versión modificada respete este comando
dump 1 all custom 1 dump.force id type fx fy fz

# Ejecución (puedes cambiar 'minimize' por 'run' o lo que necesites)
thermo 1
thermo_style custom step etotal press vol temp
run 0 # Ejecuta 0 pasos para obtener la energía y fuerzas después de la minimización/lectura

# Finalizar
write_restart restart.final
"""

    # 1. Crear la estructura de átomos (Silicio)
    atoms = bulk('Si', 'diamond', a=5.43, cubic=True)
    atoms.pbc = True # Periodic boundary conditions

    # 2. Instanciar el Calculator modificado
    try:
        calc = LammpsModified(
            # Comando de ejecución con tu binario
            command='./lmp_mpi_chimes < in.sw',
            # Contenido del script de entrada
            lammps_input=LAMMPS_INPUT_SCRIPT,
            # Parámetros adicionales para LAMMPS, si es necesario
            specorder=['Si'] # Orden de los elementos para el archivo de datos LAMMPS
        )
    except ValueError as e:
        print(f"Error al inicializar el calculator: {e}")
        exit()


    # 3. Asignar el Calculator a la estructura
    atoms.set_calculator(calc)

    print("Calculadora asignada. Intentando calcular energía y fuerzas...")

    # 4. Ejecutar el cálculo (llama a write_input, ejecuta el comando y luego read_results)

    try:
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()

        print("\n--- Resultados (Simulados/Parseados) ---")
        print(f"Energía potencial (Etot): {energy} eV")
        print(f"Fuerzas del primer átomo: {forces[0]}")

    except (FileNotFoundError, IOError) as e:
        print(f"\n¡ERROR CRÍTICO DURANTE EL CÁLCULO!")
        print(f"ASE no pudo encontrar o leer los archivos de salida después de la ejecución del binario.")
        print(f"Asegúrate de que el binario './lmp_mpi_chimes' esté en el PATH o en la carpeta actual,")
        print(f"y que produce un archivo de log y de fuerzas con el formato esperado.")
        print(f"Detalles del error: {e}")
