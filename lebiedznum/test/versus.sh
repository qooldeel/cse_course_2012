#!/bin/bash 

FILE=$1
DIM=$2
NAME=$3

#Y1=0.000032
#Y2=0.000037

LW=2

gnuplot -persist << PLOT
set xlabel 't'
set ylabel 'y(t)'
set title "${NAME}"
#set yrange [${Y1}:${Y2}]
plot for [i=2:${DIM}+1] '${FILE}' using 1:i with lines lw ${LW} notitle
#plot '${FILE}' using 1:3 with lines lw ${LW} title 'y_2'

quit
PLOT