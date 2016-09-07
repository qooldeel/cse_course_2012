#!/bin/bash


gnuplot -persist << PLOT
a = 1
b = 3

set terminal epslatex size 8.5cm,7.5cm color colortext
set output "trappingreg.tex"
set format xy "$%g$"

yfix = b/a
set title "Hauptisoklinen des Brusselators mit \$a, b = \\\\mathrm{const.}\$, Trapping-Region (Poincar\\\'e-Bendixson)"
set label 1 " " at 1.,yfix point pt 7 ps 1 
set label 2 "fixed point" at 0.9,b+0.75

set xlabel "\$x\$"
set ylabel "\$y\$"

set arrow from 0.4,7 to 0.5,7  ## horizontal on y'=0
set arrow from 0.8,7 to 1.3,4   ## x' >0, y' < 0
set arrow from 1.3,3.5 to 1.3,2 ## vertical on x' = 0
set arrow from 1.4,2.23 to 1.3,1.5   ## x' < 0, y' < 0
set arrow from 1.22,2.65 to 1,2.65  ## horizontal on y'=0
set arrow from 0.7,2.65 to 0.6,3 ## x' < 0, y' > 0
set arrow from 0.4,3 to 0.4,6  ## vertical on x' = 0

set xrange[0:1.5]
set yrange[0:9]
plot ((b+1)*x-1)/(a*x**2) lw 2 title "$\\\\dot{x} = 0$", (b*x)/(a*x**2) lw 2 title "$\\\\dot{y} = 0$"
quit
PLOT