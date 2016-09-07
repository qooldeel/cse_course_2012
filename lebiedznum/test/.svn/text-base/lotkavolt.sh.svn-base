#!/bin/bash 

FILE=$1

XLAB=Hasen
YLAB=Füchse

gnuplot -persist << PLOT
set xlabel '${XLAB}'
set ylabel '${YLAB}'
plot '${FILE}' using 2:3 with lines lw 2 title "Lotka-Volterra"

quit
PLOT