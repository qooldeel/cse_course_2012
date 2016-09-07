#!/bin/bash 

FILE=$1
MU=$2

gnuplot -persist << PLOT
set xlabel 'y1'
set ylabel 'y2'
set title "Van der Pol Oscillator: mu = ${MU}"
plot '${FILE}' using 2:3 with lines lw 2 title "phase portrait"


quit
PLOT