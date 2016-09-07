%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% cycle.m
% cycle graphs: main
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-12-07
% last change: 2012-12-07
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function cycle

% Anfangswerte
r=5.01;
theta=270;

%Zeit
tf=10000;

% Integriere das Modell mit Anfangswerten und Parametern    
[t_out,x_out] = integrate(r,theta,tf);

% Rufe die Plotfunktion
%plot_versus_time(t_out,x_out);

% Phasenportrait
plot_portrait(t_out,x_out);

return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
