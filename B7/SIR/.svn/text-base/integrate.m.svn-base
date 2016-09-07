%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% integrate.m
% SIR-Modell: Integratoraufruf
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-11-26
% last change: 2012-11-27
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function [t_out,x_out]=integrate(S,I,R,beta,r,tf)

% Integratoroptionen
options = odeset('RelTol', 1e-8,'AbsTol',1e-10,'Refine',1);

% Zeitintervall
t0 = 0;

x( 1) = S;
x( 2) = I;
x( 3) = R;

% Rufe Integrator
[t_out,x_out] = ode45(@(t,x)rhs(t,x,beta,r),[t0,tf],x,options);

return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
