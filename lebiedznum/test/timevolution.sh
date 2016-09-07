#!/bin/bash

FILE=$1

gnuplot -persist << PLOT
set xlabel 'time t'
set ylabel 'y(t)'
plot '${FILE}' using 1:2 with lines lw 2 title 'y1(t)', '${FILE}' using 1:3 with lines lw 2 title 'y2(t)'

quit
PLOT