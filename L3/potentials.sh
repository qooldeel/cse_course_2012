#!/bin/bash 

gnuplot -persist << PLOT

set xzeroaxis
set yzeroaxis

set term x11 0    # open 1st window
f(x) = 1./3.*x**3 - 0.5*x**2
set xrange[-1.75:1.75]
plot f(x) lw 2 title 'V(x)=1./3x**3 - 0.5x**2 von Aufg. 3.1 a)' 

set term x11 1 #open 2nd window
set xrange [-2.:2]
plot cosh(x)  lw 2 title 'V(x)=cosh(x) von Aufg. 3.1 b)'

set term x11 2 #open 3rd window
f(x) = 0.25*x**4 - 0.5*x**2
set xrange[-1.5:1.5]
plot f(x) lw 2 title 'V(x)= 0.24x**4-0.5*x**2 von Aufg. 3.1 c)' 

quit
PLOT