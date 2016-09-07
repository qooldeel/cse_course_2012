#!/bin/bash 

FILE=$1

XLAB=Hasen
YLAB=Füchse

gnuplot -persist << PLOT
set xlabel '${XLAB}'
set ylabel '${YLAB}'
plot '${FILE}' using 2:3:4:5 with vectors head filled lt 2 title "direction field"

quit
PLOT