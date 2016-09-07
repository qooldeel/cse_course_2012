#!/bin/bash 

FILE=$1
MU=$2

gnuplot -persist << PLOT
set xlabel 'x'
set ylabel 'y'
set title "Easy Hopf bifurcation example: mu = ${MU}"
plot '${FILE}' using 2:3 with lines lw 2 title "phase portrait"


quit
PLOT