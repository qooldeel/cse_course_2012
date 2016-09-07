%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% rhs.m
% SIR-Modell: echte Seite
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-11-26
% last change: 2012-11-27
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function xdot=rhs(t,x,beta,r)

% x(1)  S
% x(2)  I
% x(3)  R

% Vektor anlegen für Geschwindigkeit
xdot(1:3,1)=0.0;

% Model
xdot(1,1) = - beta * x(1) * x(2);
xdot(2,1) =   beta * x(1) * x(2) - r * x(2);
xdot(3,1) =                        r * x(2);

% Verlust der Immunität
a = 0.25;
xdot(1,1) = xdot(1,1) + a*x(3);
xdot(3,1) = xdot(3,1) - a*x(3);

return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
