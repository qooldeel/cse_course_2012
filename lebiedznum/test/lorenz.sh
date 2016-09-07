#!/bin/bash 

FILE=$1

LW=0.5

gnuplot -persist << PLOT
set xlabel 'x'
set ylabel 'y'
set zlabel 'z'
set title "Lorenz attractor"
splot '${FILE}' using 2:3:4 with lines lw ${LW} title "phase portrait"

quit
PLOT