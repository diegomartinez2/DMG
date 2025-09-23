#!/bin/bash
analyze_phonons.py --calc tau --temp 300 si222.result > tau300K_10.dat
echo "$ gnuplot
gnuplot> set xrange [1:]
gnuplot> set logscale y
gnuplot> set xlabel "Phonon frequency (cm^{-1})"
gnuplot> set ylabel "Phonon lifetime (ps)"
gnuplot> plot "tau300K_10.dat" using 3:4 w p"
