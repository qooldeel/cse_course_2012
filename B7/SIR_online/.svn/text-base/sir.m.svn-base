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
% daraus ergibt sich für den Rest



% Impfung mit Impfrate v (H)
%...

%Zeit tf
%...

% Parameter
%...


% Integriere das Modell mit Anfangswerten und Parametern    
[t_out,x_out] = integrate(S,I,R,beta,r,tf);

% Rufe die Plotfunktion
plot_people(t_out,x_out);

% Phasenportrait
plot_portrait(N,beta,r,t_out,x_out);

return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
