#!/bin/bash

FILE=$1
A=$2
B=$3

gnuplot -persist << PLOT
set title "Brusselator"
set xrange[0:4]
set yrange[0:7]
plot '${FILE}' using 2:3 with lines lw 2 title "phase portrait", '${FILE}' using 2:4 with lines lw 2 title "x' = 0", '${FILE}' using 2:5 with lines lw 2 title "y' = 0"

quit
PLOT