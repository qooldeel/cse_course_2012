%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% plot_versus_time.m
% cycle Modell: plottet Werte gegen die Zeit
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-12-07
% last change: 2012-12-07
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function plot_versus_time(t_out,x_out)

%Plotte gegen Zeit
figure; hold;

plot(...
    t_out(:,1),x_out(:,1), ...
    t_out(:,1),x_out(:,2) ...
    );

leg = legend('r','theta');
xlabel('time');
ylabel('variables');
  
return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
