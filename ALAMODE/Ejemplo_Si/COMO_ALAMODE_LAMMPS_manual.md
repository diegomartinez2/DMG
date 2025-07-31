Es crucial entender que este formato es **característico de versiones más antiguas o compilaciones específicas de ALAMODE** (o su precursor, el módulo `alm` de Phonopy) y difiere significativamente de la sintaxis actual.

-----

# Manual de Usuario: Cálculo de Fonones con ALAMODE (Versión Antigua) y LAMMPS

Este manual describe el proceso paso a paso para calcular las propiedades fonónicas del Silicio utilizando tu versión de ALAMODE para el ajuste de constantes de fuerza y los cálculos de fonones, y LAMMPS para la obtención de las fuerzas atómicas.

**Basado en el formato de entrada de ALAMODE proporcionado:**

  * Los bloques se cierran y el siguiente se abre con `/&NOMBRE_DEL_BLOQUE`.
  * `NAT`, `NKD`, `KD`, `NORDER` se definen en el bloque `&general`.
  * La celda unitaria se define con un factor de escala global (en Bohr) y vectores de la red en `&cell`.
  * Las posiciones de todos los átomos (de la supercelda) se definen en un bloque `&position` separado.
  * Las distancias de corte se definen en un bloque `&cutoff` separado.

**Requisitos Previos:**

  * **ALAMODE:** Tu compilación específica de ALAMODE.
  * **LAMMPS:** Compilado e instalado con el paquete `USER-PHONON` habilitado.
  * **Python:** Con la librería `ase` (Atomic Simulation Environment) instalada (`pip install ase`) para la conversión de formatos y construcción de superceldas. También para `alm_extract_lammps.py`.
  * **Potencial Interatómico:** Un archivo de potencial interatómico para Silicio compatible con LAMMPS (ej. `Si.tersoff`, `Si.eam`).
  * **`phonopy` (opcional pero recomendado):** Para generar la supercelda inicial de manera sencilla. (`pip install phonopy`)

-----

### Paso 0: Preparar la Celda Unidad y Generar la Supercelda

Tu ALAMODE requiere que la estructura completa de la supercelda (con todos sus átomos) esté incrustada en el archivo de entrada. Para el Silicio, usaremos una supercelda $2 \\times 2 \\times 2$ de la celda convencional de 8 átomos, lo que resulta en una supercelda de 64 átomos.

**0.1. Archivo de Celda Unidad `Si.POSCAR` (Referencia)**
Este es el punto de partida (celda convencional de 8 átomos).

```
Silicon (Si) Diamond Cubic - Conventional Cell (a=5.43 Angstroms)
   5.4300000000000000
    1.0000000000000000    0.0000000000000000    0.0000000000000000
    0.0000000000000000    1.0000000000000000    0.0000000000000000
    0.0000000000000000    0.0000000000000000    1.0000000000000000
Si
   8
Direct
 0.0000000000000000  0.0000000000000000  0.0000000000000000
 0.2500000000000000  0.2500000000000000  0.2500000000000000
 0.0000000000000000  0.5000000000000000  0.5000000000000000
 0.2500000000000000  0.7500000000000000  0.7500000000000000
 0.5000000000000000  0.0000000000000000  0.5000000000000000
 0.7500000000000000  0.2500000000000000  0.7500000000000000
 0.5000000000000000  0.5000000000000000  0.0000000000000000
 0.7500000000000000  0.7500000000000000  0.2500000000000000
```

**0.2. Generar la Supercelda de 64 Átomos (`SPOSCAR`)**
Usaremos `phonopy` para esto, es la forma más sencilla.

1.  Asegúrate de que el archivo de tu celda unidad se llame `POSCAR` (tal cual, sin `Si.`).
2.  Crea un archivo llamado `phonopy_disp.conf` con el siguiente contenido:
    ```
    DIM = 2 2 2
    ```
3.  Ejecuta `phonopy` en tu terminal:
    ```bash
    phonopy -c POSCAR phonopy_disp.conf
    ```
    Esto generará un archivo llamado `SPOSCAR`, que contiene la descripción de tu supercelda $2 \\times 2 \\times 2$ (64 átomos). Este `SPOSCAR` es el que usaremos para copiar las coordenadas a los archivos de entrada de ALAMODE.

**0.3. Determinar el Factor de Escala de la Supercelda en Bohr:**
La constante de red de la celda convencional de Si es $5.43$ Å. Para una supercelda $2 \\times 2 \\times 2$, la constante de red de la supercelda será $2 \\times 5.43 = 10.86$ Å.
Tu ALAMODE requiere este valor en **Bohr**.
Factor de conversión: $1 \\text{ Bohr} \\approx 0.529177 \\text{ Å}$.
Constante de red de la supercelda en Bohr: $10.86 \\text{ Å} / 0.529177 \\text{ Å/Bohr} \\approx 20.5212 \\text{ Bohr}$.
Usaremos `20.5212` como el factor de escala en el bloque `&cell`.

-----

### Paso 1: Generar Patrones de Desplazamiento con ALAMODE (`alm_suggest`)

Crea un archivo llamado `alm_suggest.in` con el siguiente contenido. Asegúrate de **copiar las 64 líneas de coordenadas atómicas del `SPOSCAR`** generado en el Paso 0.2 y pegarlas en el bloque `&position`. El tipo de átomo para Silicio es `1` en este formato.

**`alm_suggest.in`**

```ini
! alm_suggest.in (Formato adaptado a tu manual)
&general
  PREFIX = Si_displacements_patterns ! Prefijo para los archivos de salida
  MODE = suggest              ! Modo para generar patrones de desplazamiento
  NAT = 64                    ! Número total de átomos en la supercelda
  NKD = 1                     ! Número de tipos de átomos
  KD = Si
  /

  &interaction       ! Nombre del tipo de átomo, seguido del delimitador del siguiente bloque

  NORDER = 2  # 1: armónico, 2: cúbico  ! Orden máximo de las FC a considerar
  /
  &cell
!--------------------------------------------------------------------------------------------------
20.5212                     ! Factor de escala de la red en Bohr (Constante de red de la supercelda 2x2x2)
1.0 0.0 0.0                 ! Vectores de la red normalizados
0.0 1.0 0.0
0.0 0.0 1.0
/

&cutoff
!--------------------------------------------------------------------------------------------------
Si-Si 5.0 3.0                  ! Distancia de corte para interacciones armónicas (FC2) en Angstroms y Distancia de corte para interacciones cúbicas (FC3) en Angstroms
# *-* 5.0 3.0                ! También se pueden poner asteriscos
/

&position
!--------------------------------------------------------------------------------------------------
! ATENCIÓN: PEGA AQUÍ LAS 64 LÍNEAS DE POSICIONES ATÓMICAS DE TU ARCHIVO SPOSCAR.
! El formato es: [ÍNDICE_TIPO_ÁTOMO] [COORD_X] [COORD_Y] [COORD_Z]
! Para Si, el índice de tipo de átomo es '1'. Las coordenadas son fraccionarias.
! Ejemplo de las primeras 8 (las 64 deben ir aquí):
1 0.0000000000000000 0.0000000000000000 0.0000000000000000
1 0.1250000000000000 0.1250000000000000 0.1250000000000000
1 0.0000000000000000 0.2500000000000000 0.2500000000000000
1 0.1250000000000000 0.3750000000000000 0.3750000000000000
1 0.2500000000000000 0.0000000000000000 0.2500000000000000
1 0.3750000000000000 0.1250000000000000 0.3750000000000000
1 0.2500000000000000 0.2500000000000000 0.0000000000000000
1 0.3750000000000000 0.3750000000000000 0.1250000000000000
! ... (CONTINÚA CON LAS 56 LÍNEAS RESTANTES DE TU SPOSCAR)
! Asegúrate de que no haya un delimitador /& al final del archivo si no hay más bloques.
```

**Ejecución:**

```bash
alm alm_suggest.in
```

**Archivos de Salida Clave:**

  * `Si_displacements_patterns.pattern`: Contiene los patrones de desplazamiento que debes simular en LAMMPS.
  * `Si_displacements_patterns.xyz`: Una representación visual de la supercelda original (sin desplazamientos).

-----

### Paso 2: Ejecutar Simulaciones con LAMMPS para Obtener las Fuerzas

**2.1. Crear un Script de LAMMPS Genérico (`in.lammps_template`)**

Este script será modificado y ejecutado para cada patrón de desplazamiento.

```lammps
# in.lammps_template - LAMMPS input script for ALAMODE force calculations

# Initialization
units metal
atom_style atomic
boundary p p p

# Crea los átomos de la supercelda a partir de un archivo de datos LAMMPS
read_data Si_supercell.data

# Potencial (Ejemplo: Tersoff para Si)
pair_style tersoff
pair_coeff * * Si.tersoff Si

# Configuraciones de vecinos
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes

# Fija los átomos en sus posiciones desplazadas y calcula las fuerzas
# El fix phonon usa el archivo .pattern de ALAMODE para aplicar el desplazamiento.
fix 1 all phonon 1 0 Si_displacements_patterns.pattern PLACEHOLDER_PATTERN_INDEX

# Salida de fuerzas y posiciones para ALAMODE
dump 1 all custom 1 dump.all id type x y z fx fy fz
dump_modify 1 sort id

# Ejecuta un solo paso de tiempo para calcular las fuerzas
run 0
```

**2.2. Script Bash para Generar `Si_supercell.data` y Ejecutar LAMMPS**

Este script convertirá el `SPOSCAR` a un formato `.data` de LAMMPS y luego ejecutará LAMMPS para cada patrón.

```bash
#!/bin/bash

# --- Configuración ---
LAMMPS_EXE="lmp_mpi" # O lmp_serial, o la ruta completa a tu binario LAMMPS
N_CORES=72           # Ajusta esto a los núcleos que quieras usar para MPI
PATTERN_FILE="Si_displacements_patterns.pattern"
TEMPLATE_LAMMPS_IN="in.lammps_template"
POTENTIAL_FILE="Si.tersoff" # Asegúrate de que este archivo exista

# --- Limpieza de archivos anteriores ---
rm -f dump.* forces_lammps_run_*.dat log.lammps_* Si_supercell.data

echo "--- Paso 2.1: Convertir SPOSCAR (supercelda) a LAMMPS data file ---"
# Utiliza ASE para convertir el SPOSCAR (64 átomos) a un archivo .data para LAMMPS
# Asegúrate de que 'ase' esté instalado: pip install ase
python -c "from ase.io import read, write; atoms = read('SPOSCAR'); write('Si_supercell.data', atoms, format='lammps-data')"
if [ $? -ne 0 ]; then
    echo "ERROR: Falló la conversión de SPOSCAR a Si_supercell.data. Asegúrate de que ASE esté instalado y SPOSCAR sea válido."
    exit 1
fi
echo "SPOSCAR (supercelda) convertido a Si_supercell.data"

# Copia el potencial al directorio de trabajo
cp $POTENTIAL_FILE .

echo "--- Paso 2.2: Generando estructuras y ejecutando LAMMPS para cada patrón de desplazamiento ---"

NUM_PATTERNS=$(grep -c "BEGIN_PATTERN" ${PATTERN_FILE})
echo "Número total de patrones de desplazamiento a simular: ${NUM_PATTERNS}"

for i in $(seq 0 $((NUM_PATTERNS-1)))
do
    echo "Simulando patrón ${i}..."

    LAMMPS_IN_FILE="in.lammps_run_${i}"

    # Crea un script LAMMPS temporal para cada corrida
    cp "${TEMPLATE_LAMMPS_IN}" "${LAMMPS_IN_FILE}"

    # Modifica el fix phonon para este patrón específico
    sed -i "s/PLACEHOLDER_PATTERN_INDEX/${i}/" "${LAMMPS_IN_FILE}"

    # Renombra el archivo de dump para cada corrida
    sed -i "s/dump.all/dump.all_${i}/" "${LAMMPS_IN_FILE}"

    # Ejecuta LAMMPS (ajusta 'mpirun' si usas MPI, o si es lmp_serial)
    echo "Running LAMMPS for pattern ${i} with command: mpirun -np ${N_CORES} ${LAMMPS_EXE} -in ${LAMMPS_IN_FILE}"
    mpirun -np ${N_CORES} ${LAMMPS_EXE} -in "${LAMMPS_IN_FILE}" > "log.lammps_${i}"

done

echo "Todas las simulaciones LAMMPS completadas."

echo "--- Paso 2.3: Extraer Fuerzas de LAMMPS con alm_extract_lammps.py ---"
# alm_extract_lammps.py procesará todos los archivos dump.all_X y el archivo .pattern
# para generar un único archivo de fuerzas en el formato que alm_optimize espera.
alm_extract_lammps.py --prefix Si_displacements_patterns \
                      --lammps_dump "dump.all_*" \
                      --output Si_LAMMPS_extracted_forces.pattern

echo "Fuerzas extraídas y guardadas en Si_LAMMPS_extracted_forces.pattern"

# Limpia los archivos dump.all_* si no los necesitas más
# rm -f dump.all_*
```

-----

### Paso 3: Optimizar Constantes de Fuerza con ALAMODE (`alm_optimize`)

Crea un archivo llamado `alm_optimize.in`. **Copia exactamente el mismo contenido de `&cell`, `&cutoff` y `&position` que usaste en `alm_suggest.in`**.

**`alm_optimize.in`**

```ini
! alm_optimize.in (Formato adaptado a tu manual)
&general
  PREFIX = Si_optimize        ! Prefijo para los archivos de salida
  MODE = optimize             ! Modo para optimizar las constantes de fuerza
  NAT = 64                    ! Número total de átomos en la supercelda
  NKD = 1
  KD = Si
  /

  &interaction

  NORDER = 2
  /

  &cell           ! Orden de las FC a optimizar
!--------------------------------------------------------------------------------------------------
20.5212                     ! Factor de escala de la red en Bohr (COPIAR del alm_suggest.in)
1.0 0.0 0.0                 ! Vectores de la red normalizados (COPIAR del alm_suggest.in)
0.0 1.0 0.0
0.0 0.0 1.0
/

&cutoff
!--------------------------------------------------------------------------------------------------
Si-Si 5.0 3.0                  ! Distancia de corte para FC2 (COPIAR del alm_suggest.in)
# *-* 5.0 3.0                ! Distancia de corte para FC3 (COPIAR del alm_suggest.in)
/

&position
!--------------------------------------------------------------------------------------------------
! ATENCIÓN: PEGA AQUÍ LAS 64 LÍNEAS DE POSICIONES ATÓMICAS DE TU ARCHIVO SPOSCAR.
! DEBEN SER IDÉNTICAS A LAS DE alm_suggest.in
! Ejemplo de las primeras 8:
1 0.0000000000000000 0.0000000000000000 0.0000000000000000
1 0.1250000000000000 0.1250000000000000 0.1250000000000000
! ... (62 LÍNEAS RESTANTES)
!--------------------------------------------------------------------------------------------------
&optimize                     ! Parámetros de optimización
  DISPLACEMENT_IMAGES = Si_LAMMPS_extracted_forces.pattern ! Archivo con los desplazamientos y fuerzas
  ALGORITHM = iterative
  LINEAR_MODEL = 1
  FORCE_CONV_THRESHOLD = 1.0e-4
  MAX_ITERATION = 100
/
```

**Ejecución:**

```bash
alm alm_optimize.in
```

**Archivos de Salida Clave:**

  * `Si_LAMMPS.FC2.fcs`: Constantes de fuerza armónicas (2º orden).
  * `Si_LAMMPS.FC3.fcs`: Constantes de fuerza cúbicas (3º orden).
  * `Si_LAMMPS.alm`: Archivo de log detallado del proceso de optimización.

-----

### Paso 4: Calcular Propiedades Fonónicas con ALAMODE (`anphon`)

Crea un archivo llamado `anphon.in`. **Copia de nuevo el mismo contenido de `&cell`, `&cutoff` y `&position`**.

**`anphon.in`**

```ini
! anphon.in (Formato adaptado a tu manual)
&general
  PREFIX = Si_LAMMPS_Phonons   ! Prefijo para los archivos de salida
  MODE = phonon               ! Modo para calcular propiedades fonónicas
  NAT = 64                    ! Número total de átomos en la supercelda
  NKD = 1
  KD = Si
  /

  &interaction

  NORDER = 2
  /

  &cell           ! Orden de las FC utilizadas (armónicas y cúbicas)
!--------------------------------------------------------------------------------------------------
20.5212                     ! Factor de escala de la red en Bohr (COPIAR del alm_suggest.in)
1.0 0.0 0.0                 ! Vectores de la red normalizados (COPIAR del alm_suggest.in)
0.0 1.0 0.0
0.0 0.0 1.0
/

&cutoff
!--------------------------------------------------------------------------------------------------
Si-Si 5.0 3.0                  ! Cutoff used when FCs were generated (COPIAR del alm_suggest.in)
# *-* 5.0 3.0                ! Cutoff used when FCs were generated (COPIAR del alm_suggest.in)
/

&position
!--------------------------------------------------------------------------------------------------
! ATENCIÓN: PEGA AQUÍ LAS 64 LÍNEAS DE POSICIONES ATÓMICAS DE TU ARCHIVO SPOSCAR.
! DEBEN SER IDÉNTICAS A LAS DE alm_suggest.in
! Ejemplo de las primeras 8:
1 0.0000000000000000 0.0000000000000000 0.0000000000000000
1 0.1250000000000000 0.1250000000000000 0.1250000000000000
! ... (62 LÍNEAS RESTANTES)
!--------------------------------------------------------------------------------------------------
&phonon                     ! Parámetros de cálculo de fonones
  FC_FILE = Si_LAMMPS.FC2.fcs ! Ruta al archivo de constantes de fuerza armónicas
  FC_FILE = Si_LAMMPS.FC3.fcs ! Ruta al archivo de constantes de fuerza cúbicas

  NONANALYTIC = 1             ! Corrección no analítica (generalmente se incluye)
  DOS_BIN_WIDTH = 0.05        ! Ancho de la celda (bin) para DOS en THz
  TEMP_MIN = 0.0              ! Temperatura mínima (K)
  TEMP_MAX = 300.0            ! Temperatura máxima (K)
  TEMP_INTERVAL = 25.0        ! Intervalo de temperatura (K)

  ! SUPERCELL: Este parámetro en anphon.in *no* define la estructura de entrada,
  ! sino el tamaño de la supercelda para la diagonalización de la matriz dinámica.
  ! Podría ser redundante si la estructura ya está en la supercelda.
  ! Si da error, puedes probar a comentarlo o quitarlo. Lo común es que vaya aquí.
  SUPERCELL = 4 4 4           ! Ejemplo: Supercelda 4x4x4 para el cálculo de la matriz dinámica.

  ISMEAR = 0                  ! Método de ensanchamiento (Gaussiano)
  SMEAR = 0.005               ! Parámetro de ensanchamiento (Ry)
/

&kpoint                     ! Parámetros de puntos K para la estructura de banda y DOS
  KPOINT_MODE = BAND
  KPOINT_MESH = 16 16 16    ! Malla de k-points para DOS y conductividad térmica (para kappa, necesitarás una malla más densa)

  ! Ruta de k-points para la dispersión de bandas (para Silicio tipo FCC)
  G 0.000000 0.000000 0.000000 X 0.500000 0.500000 0.000000 51
  X 0.500000 0.500000 0.000000 W 0.500000 0.250000 0.000000 51
  W 0.500000 0.250000 0.000000 L 0.500000 0.500000 0.500000 51
  L 0.500000 0.500000 0.500000 K 0.375000 0.750000 0.375000 51
  K 0.375000 0.750000 0.375000 G 0.000000 0.000000 0.000000 51
/

&dos
  DOS_CALC = 1
/

&thermal
  THERMAL_COND = 1
  NKS_THERMAL = 1000          ! Puntos aleatorios para conductividad. Para precisión, ¡mucho mayor!
/
```

-----

¡Este manual debería ser mucho más preciso para tu entorno ALAMODE\! Alinea cada parámetro con el formato que tu binario espera. Recuerda que la clave es la consistencia en la definición de la supercelda de 64 átomos en todos los archivos de entrada de ALAMODE.
