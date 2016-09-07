%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% integrate.m
% cycle-Modell: Integratoraufruf
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-12-07
% last change: 2012-12-07
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function [t_out,x_out]=integrate(r,theta,tf)

% Integratoroptionen
options = odeset('RelTol', 1e-12,'AbsTol',1e-14,'Refine',1);

% Zeitintervall
t0 = 0;

x( 1) = r;
x( 2) = theta;

% Rufe Integrator
[t_out,x_out] = ode45(@(t,x)rhs(t,x),[t0,tf],x,options);

return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
