echo "1. Otro análisis para la conductividad termica acumulada:"
analyze_phonons.py --calc cumulative --temp 300 --length 10000:5 si222.result > cumulative_300K_10.dat
echo "2. Dibujamos:"
echo "$ gnuplot"
echo "gnuplot> set logscale x"
echo "gnuplot> set xlabel "L (nm)""
echo "gnuplot> set ylabel "Cumulative kappa (W/mK)""
echo "gnuplot> plot "cumulative_300K_10.dat" using 1:2 w lp""
