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
