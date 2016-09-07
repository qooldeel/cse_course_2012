#!/bin/bash 

FILE=$1
TITLE=$2
NUMSOL=3
EXACTSOL=5

gnuplot -persist << PLOT
set xlabel 't'
set ylabel 'y2(t)'
set title "${TITLE}"
plot '${FILE}' using 1:${NUMSOL} with lines lw 2 title "numerical solution", '${FILE}' using 1:${EXACTSOL} with lines lw 2 title "exact solution"

quit
PLOT