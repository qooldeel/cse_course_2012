%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% integrate.m
% FHN-Modell: Integratoraufruf
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-12-03
% last change: 2012-12-03
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function [t_out,x_out]=integrate(u,v,epsilon,gamma,delta,a,I,tf)

% Integratoroptionen
options = odeset('RelTol', 1e-12,'AbsTol',1e-14,'Refine',1);

% Zeitintervall
t0 = 0;

x( 1) = u;
x( 2) = v;

% Rufe Integrator
[t_out,x_out] = ode45(@(t,x)rhs(t,x,epsilon,gamma,delta,a,I),[t0,tf],x,options);

return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
