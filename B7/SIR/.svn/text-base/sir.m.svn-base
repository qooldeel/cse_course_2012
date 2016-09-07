%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% sir.m
% SIR-Modell: main
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-11-26
% last change: 2012-11-27
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function sir

% Anfangswerte: Einer von 500 ist infiziert
N = 500;
I = 1;

% daraus ergibt sich für den Rest
S = N-I;
R = N-S-I;

% Impfung mit Impfrate v
v = 0.0;
S = S - v*N;
R = R + v*N;

%Zeit
tf=100;

% Parameter
r    = 0.1;
beta = 0.001;

% Integriere das Modell mit Anfangswerten und Parametern    
[t_out,x_out] = integrate(S,I,R,beta,r,tf);

% Rufe die Plotfunktion
plot_people(t_out,x_out);

% Phasenportrait
plot_portrait(N,beta,r,t_out,x_out);

return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
