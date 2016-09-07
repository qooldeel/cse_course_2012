%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% plot_portrait.m
% FHN-Modell: plottet Phasenportrait
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-12-03
% last change: 2012-12-03
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function plot_portrait(epsilon,gamma,delta,a,I,t_out,x_out)

% Gitter
u_min = -1.5;
v_min = -1.5;
u_max = 1.5;
v_max = 1.5;

[U,V] = meshgrid(u_min:0.1:u_max,v_min:0.1:v_max); 


% Berechnung der rechten Seite
t=0;

for j=1:size(U,1)
    for k=1:size(U,2)
        x(1) = U(j,k);
        x(2) = V(j,k);
        xdot = rhs(t,x,epsilon,gamma,delta,a,I);
       
        Udot(j,k) = xdot(1,1);
        Vdot(j,k) = xdot(2,1);
    end;
end;

% quiverplot
figure; hold;
xlim([u_min,u_max]);
ylim([v_min,v_max]);
xlabel('u');
ylabel('v');
quiver(U,V,Udot,Vdot);

% und die Nullklinen
x=linspace(u_min,u_max,1000);
y(1:length(x))=I + x(1:length(x)).*(a-x(1:length(x))).*(x(1:length(x))-1);
plot(x,y,'color','green');
y(1:length(x))=(x(1:length(x)) + delta)/gamma;
plot(x,y,'color','red');

% und noch die Trajektorie
plot( x_out(:,1), x_out(:,2), 'color', 'black' );

return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
