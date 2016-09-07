#!/bin/bash 

FILE=$1

LW=2

gnuplot -persist << PLOT
set xlabel 't'
set ylabel 'y'
set title "Test"
plot '${FILE}' using 1:2 with lines lw ${LW} title "x1", '${FILE}' using 1:3 with lines lw ${LW} title "x2"

quit
PLOT