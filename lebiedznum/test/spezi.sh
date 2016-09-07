#!/bin/bash

FILE=$1

gnuplot -persist << PLOT
##{/Symbol m}x_i is only working when using set terminal postscripts eps enhanced 
set xlabel "mu**2 - 4"
set ylabel "mu"
set xzeroaxis
set term x11 0
plot '${FILE}' using 2:1 with lines lw 2 notitle
set term x11 1
set yzeroaxis
set xlabel "real: a"
set ylabel 'imag: b'
set title "eigenvalues lambda = a + bi"
plot '${FILE}' using 3:4 with linespoints lw 2 pt 4 title "lambda_1", '${FILE}' using 5:6 with linespoints lw 2 pt 5 title "lambda_2"
quit
PLOT