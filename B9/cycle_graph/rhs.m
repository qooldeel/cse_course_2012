%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% rhs.m
% cycle Modell: rechte Seite
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-12-07
% last change: 2012-12-07
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function xdot=rhs(t,x,epsilon,gamma,delta,a,I)

% x(1)  r
% x(2)  theta

% Vektor anlegen für Geschwindigkeit
xdot(1:2,1)=0.0;

% Model
xdot(2,1) = x(1)^2*(sin(x(2)))^2 + (x(1)^2*(cos(x(2)))^2 - 1.0)^2 ;
xdot(1,1) = x(1) * (1.0 - x(1)^2) * xdot(2,1);

return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
