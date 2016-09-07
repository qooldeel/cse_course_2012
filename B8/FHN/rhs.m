%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% rhs.m
% FHN-Modell: echte Seite
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-12-03
% last change: 2012-12-03
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function xdot=rhs(t,x,epsilon,gamma,delta,a,I)

% x(1)  u
% x(2)  v

% Vektor anlegen für Geschwindigkeit
xdot(1:2,1)=0.0;

% Model
f = x(1) * (a - x(1)) * (x(1) - 1.0);

xdot(1,1) = f - x(2) + I;
xdot(2,1) = epsilon*(x(1) - gamma*x(2) + delta);

return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
