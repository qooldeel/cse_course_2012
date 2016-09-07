%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% fhn.m
% FHN-Modell: main
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-12-03
% last change: 2012-12-03
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function fhn

% Anfangswerte
u=-0.682;
v=-0.365;

%Zeit
tf=500;

% Parameter
epsilon=0.01;
gamma=0.5;
delta=0.0;
delta=0.5;
a=-1.0;

% Strom
I=1.02;

% Integriere das Modell mit Anfangswerten und Parametern    
[t_out,x_out] = integrate(u,v,epsilon,gamma,delta,a,I,tf);

% Rufe die Plotfunktion
plot_versus_time(t_out,x_out);

% Phasenportrait
plot_portrait(epsilon,gamma,delta,a,I,t_out,x_out);

return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
