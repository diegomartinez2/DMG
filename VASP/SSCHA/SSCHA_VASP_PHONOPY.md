0) File "POSCAR":

SrTiO3 Unit cells
3.9048       <--cell parameter
Sr 0.0 0.0 0.0
O 0.0 0.5 0.5
O 0.5 0.0 0.5
O 0.5 0.5 0.0
Ti 0.5 0.5 0.5

1) phonopy -d --dim='3 3 3'   <-- create supercell by displacements
2) cp POSCAR POSCAR_UNITCELL
3) cp SPOSCAR POSCAR
4) mkdir harmonic_calculation
4.1) cp INCAR.harmonic ./harmonic/INCAR
4.2) VASP calculation
4.3) Copy results into sscha directory and go there.
--------
5) phonopy --fc vasprun.xml  <-- It will fail if the xml is incomplete.
5.1) cp POSCAR_UNITCELL POSCAR
5.2) run "interface.py":
* PHONOPHY supercell = 3  <-- For a 3x3x3 supercell
* interpolate = y
* QxQxQ = 3 <-- if we want 3x3x3 other for other interpolation.
* Ensembles = 300 <-- Number of ensembles to create. For Sobol use n**2.
* Pop = 1
* nquirr = 4 <-- The number of irreducible elements (or just look the number of dyn* files).
* TEmp = 300 <-- Temperature in Kelvin
* Positivate = n <-- or 'y' to make positive definite.
6) cp INCAR.sc ./pop1/vasp/INCAR
6.1) cp POTCAR ./pop1/vasp/POTCAR
6.2) cp ML_FF ./pop1/vasp/ML_FF
6.3) Edit run.bach:
* POPULATION = 1
* IONS = 135 <-- Number of atoms in the supercell.
* NCONFSSCHA = 300 <-- Number of ensembles in the population.
----
VASP runs
----
