%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% plot_portrait.m
% SIR-Modell: plottet Phasenportrait
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-11-27
% last change: 2012-11-27
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function plot_portrait(N,beta,r,t_out,x_out);

% Gitter
s_max = 500;
i_max = 300;

[s,i] = meshgrid(0:5:s_max,0:3:i_max); 


% Berechnung der rechten Seite
t=0;

for j=1:size(s,1)
    for k=1:size(s,2)
        x(1) = s(j,k);
        x(2) = i(j,k);
        x(3) = N - x(1) - x(2);
        xdot = rhs(t,x,beta,r);
       
        sdot(j,k) = xdot(1,1);
        idot(j,k) = xdot(2,1);
    end;
end;

% quiverplot
figure; hold;
xlim([0,s_max]);
ylim([0,i_max]);
xlabel('S');
ylabel('I');
quiver(s,i,sdot,idot);

% und noch die Trajektorie
plot( x_out(:,1), x_out(:,2), 'color', 'red' );

% und noch den Punkt
plot(r/beta,max(x_out(:,2)),'x','MarkerSize',10,'linewidth',3);
    
return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
