#!/bin/bash 

FILE=$1

gnuplot -persist << PLOT
set xlabel 'Re(z)'
set ylabel 'Im(z)'
set title "Aufgabe 11.1" 
plot '${FILE}' using 1:2 with lines lw 2 lc rgbcolor "blue" title "lambda_1", '${FILE}'  using 1:3 with lines lw 2 lc rgbcolor"orange " title "lambda_2"



quit
PLOT