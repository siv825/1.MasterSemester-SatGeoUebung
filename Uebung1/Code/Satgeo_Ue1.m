%% Satgeo Ü1
% Nadine & Ziqing

%% GOCE
a=6378137+225e3;
I=deg2rad(96.6);
e=0;
Omega=deg2rad(335);
w=deg2rad(273);
M=deg2rad(5);
GM=3.9865005*10^14;
[r,v] = kep2cart(I,Omega,w,M,e,a,GM);
r11=[r';v'];

% Atmosphärenreibung
% dc=[ ]
% [f_atm] = drag_force(dc, 225e3, v');
% input:
% 1. dc is the tabulated minimum and maximum density value (column 2 and column 3)
%    at a given altitude h_i (column 1)
% 2. h is the altitude of the satellite
% 3. v_sate is the velocity of the satellite (3*1 vector)
% output:
% f_atm is the atmospheric drag
f_atm=[1 2 3];
% r12=[r;v;f_atm];
r12=[r';v';f_atm'];
% Integration über die Bahnumläufe


options = odeset('RelTol',1e-15,'AbsTol',1e-15);

TC=2*pi*(sqrt((norm(r11(1:3)))^3/GM));  % aus dem Bachelor übernommen
t_1_sec=[0 5*TC];                       % aus dem Bachelor übernommen
[T1,Y1]=ode45(@odefun,t_1_sec,r12,options);

figure;
[x,y,z]=ellipsoid(0,0,0,6378137,6378137,6356752.3142);
surf(x, y, z);
axis equal;
grid on;
hold on;
title('Umlaufbahn GOCE')
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'LineWidth',2);

% numerische Integration
function dydt =odefun(t,y)
GM=3.9865005e14;
dydt=[y(4);y(5);y(6);-y(1)*GM/norm(y(1:3))^3+y(7);-y(2)*GM/norm(y(1:3))^3+y(8);-y(3)*GM/norm(y(1:3))^3+y(9);0;0;0]
% Ich musste dydt durch die 0;0;0 ergänzen, da r12 9 Elemente hat und somit
% dydt auh 9 Elemente braucht. Die 9 Elemente von r12 kommen 3x Position
% 3x Geschwindigkeit und 3x Drag_force
end