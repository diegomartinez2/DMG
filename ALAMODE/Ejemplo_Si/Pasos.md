# Cómo se hace...
./lmp_mpi_chimes -in lammps.in  #primero relajamos el sistema... esto nos dará 'relax.dat'
### cambiamos la estructura en alm_suggest.in con los datos de la estructura relajada "relax.dat"
alm alm_suggest.in > si_alm.log1
python displace.py --LAMMPS=relax.dat --mag=0.01 --prefix harm  -pf Si_displacements_patterns.pattern_HARMONIC
python displace.py --LAMMPS=relax.dat --mag=0.04 --prefix cubic  -pf Si_displacements_patterns.pattern_ANHARM3
./lmp_mpi_chimes < in.sw
bash run.bash
python extract.py --LAMMPS=relax.dat XFSET > DFSET_harmonic
python extract.py --LAMMPS=relax.dat XFSET.cubic* > DFSET_cubic
cat DFSET_harmonic DFSET_cubic > DFSET_merged
alm alm_optimize.in
anphon si_phband.in
plotband.py si222.bands
