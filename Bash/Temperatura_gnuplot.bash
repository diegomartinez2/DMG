#!/bin/bash
# temp_monitor.sh
sensors | grep 'Core' | awk '{print $3}' | sed 's/+//;s/°C//' | \
paste <(seq 1 8) - | \
gnuplot -p -e "
set title 'Temperaturas CPU - $(date)';
set yrange [0:100];
set style fill solid;
plot '-' with boxes lc rgb 'orange' title '°C';
pause 3
# set terminal png size 800,600;
# set output 'pressure_step.png';
# set title 'Step vs Pressure';
"
