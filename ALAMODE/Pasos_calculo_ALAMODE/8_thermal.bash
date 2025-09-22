#!/bin/bash
echo "first build the 'anphon_RTA.in' with 'Construir_alm_anphon_file.py'"
anphon anphon_RTA.in
echo "$ gnuplot
gnuplot> set logscale xy
gnuplot> set xlabel "Temperature (K)"
gnuplot> set ylabel "Lattice thermal conductivity (W/mK)"
gnuplot> plot "My_calculation_222.kl" usi 1:2 w lp"
analyze_phonons.py --calc kappa_boundary --size 1.0e+6 My_calculation_222.result > My_calculation_222_boundary_1mm.kl
echo "$ gnuplot
gnuplot> set logscale xy
gnuplot> set xlabel "Temperature (K)"
gnuplot> set ylabel "Lattice thermal conductivity (W/mK)"
gnuplot> plot "My_calculation_222.kl" usi 1:2 w lp, "My_calculation_222_boundary_1mm.kl" usi 1:2 w lp"
