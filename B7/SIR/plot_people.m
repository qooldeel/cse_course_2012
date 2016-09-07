%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% plot_people.m
% SIR-Modell: plottet Werte gegen die Zeit
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-11-26
% last change: 2012-11-27
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function plot_people(t_out,x_out)

%Plotte gegen Zeit
figure; hold;

plot(...
    t_out(:,1),x_out(:,1), ...
    t_out(:,1),x_out(:,2), ...
    t_out(:,1),x_out(:,3)  ...
    );

leg = legend('S','I','R');
xlabel('time');
ylabel('Zahl der Individuen');
  
return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
