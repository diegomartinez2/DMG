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
### conductividad terminca
anphon si_RTA.in

$ gnuplot
gnuplot> set logscale xy
gnuplot> set xlabel "Temperature (K)"
gnuplot> set ylabel "Lattice thermal conductivity (W/mK)"
gnuplot> plot "si222.kl" usi 1:2 w lp

analyze_phonons.py --calc kappa_boundary --size 1.0e+6 si222.result > si222_boundary_1mm.kl

$ gnuplot
gnuplot> set logscale xy
gnuplot> set xlabel "Temperature (K)"
gnuplot> set ylabel "Lattice thermal conductivity (W/mK)"
gnuplot> plot "si222.kl" usi 1:2 w lp, "si222_boundary_1mm.kl" usi 1:2 w lp

### Phonon lifetime
analyze_phonons.py --calc tau --temp 300 si222.result > tau300K_10.dat

$ gnuplot
gnuplot> set xrange [1:]
gnuplot> set logscale y
gnuplot> set xlabel "Phonon frequency (cm^{-1})"
gnuplot> set ylabel "Phonon lifetime (ps)"
gnuplot> plot "tau300K_10.dat" using 3:4 w p

### Cumulative thermal conductivity
analyze_phonons.py --calc cumulative --temp 300 --length 10000:5 si222.result > cumulative_300K_10.dat

$ gnuplot
gnuplot> set logscale x
gnuplot> set xlabel "L (nm)"
gnuplot> set ylabel "Cumulative kappa (W/mK)"
gnuplot> plot "cumulative_300K_10.dat" using 1:2 w lp

* nota: algo ocurre con el programa reintenta varias veces el analyze_phonons (ten en cuenta el python+cpp)

### Thermal conductivity spectrum
anphon si_RTA.in > si_RTA2.log

$ awk '{if ($1 == 300.0) print $0}' si222.kl_spec > si222_300K_10.kl_spec
$ gnuplot
gnuplot> set xlabel "Frequency (cm^{-1})"
gnuplot> set ylabel "Spectrum of kappa (W/mK/cm^{-1})"
gnuplot> plot "si222_300K_10.kl_spec" using 2:3 w l lt 2 lw 2
