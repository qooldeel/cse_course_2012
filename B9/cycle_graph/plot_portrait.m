%%.   ...1...   ...2...   ...3...   ...4...   ...5...   ...6...   ...7...   ...8
% plot_portrait.m
% cycle-Modell: plottet Phasenportrait
% Jochen Siehr
% Numerische Mathematik, Uni Ulm
% 2012-12-07
% last change: 2012-12-07
%- ----- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%

function plot_portrait(t_out,x_out)



%figure; hold;
%polar(x_out(:,2),x_out(:,1));

% Gitter1
u_min = 0;
v_min = 0;
u_max = 5;
v_max = 360;

[U,V] = meshgrid(u_min:0.05:u_max,v_min:1:v_max); 


% Berechnung der rechten Seite
t=0;

for j=1:size(U,1)
    for k=1:size(U,2)
        x(1) = U(j,k);
        x(2) = V(j,k);
        xdot = rhs(t,x);
       
        Udot(j,k) = xdot(1,1);
        Vdot(j,k) = xdot(2,1);
    end;
end;

X = U .* cos(V);
Y = U .* sin(V);
Xdot = Udot .* cos(Vdot);
Ydot = Udot .* sin(Vdot);


% Gitter2
u2_min = -5;
v2_min = -5;
u2_max = 5;
v2_max = 5;

[U2,V2] = meshgrid(u2_min:0.1:u2_max,v2_min:0.1:v2_max); 


% Berechnung der rechten Seite
t=0;

for j=1:size(U2,1)
    for k=1:size(U2,2)
        x(1) = sqrt(U2(j,k)^2 + V2(j,k)^2);
        x(2) = sign(V2(j,k))*acos(U2(j,k)/x(1));
        xdot = rhs(t,x);
       
        U2dot(j,k) = xdot(1,1);
        V2dot(j,k) = xdot(2,1);
    end;
end;

X2dot = U2dot .* cos(V2dot);
Y2dot = U2dot .* sin(V2dot);






% quiverplot
figure; hold;
%xlim([u_min,u_max]);
%ylim([v_min,v_max]);
xlabel('X');
ylabel('Y');
quiver(X,Y,Xdot,Ydot);
quiver(U2,V2,X2dot,Y2dot);

% und die Nullklinen
%x=linspace(u_min,u_max,1000);
%y(1:length(x))=I + x(1:length(x)).*(a-x(1:length(x))).*(x(1:length(x))-1);
%plot(x,y,'color','green');
%y(1:length(x))=(x(1:length(x)) + delta)/gamma;
%plot(x,y,'color','red');

% und noch die Trajektorie
%a = x_out(:,1) .* cos(x_out(:,2));
%b = x_out(:,1) .* sin(x_out(:,2));
polar(x_out(:,2),x_out(:,1),'r');


return;

%- -eof- ----- ----- ----- ----- ----- -- ----- ----- ----- ----- ----- ----- -%
