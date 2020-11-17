%% Satgeo Ü1
% Nadine & Ziqing


%% GOCE
a = 6378137+225e3; % meter
I = deg2rad(96.6); % radiant
e = 0; % keine Einheit
Omega = deg2rad(335);% radiant
w = deg2rad(273); % radiant
M = deg2rad(5); % radiant
GM = 3.9865005*10^14; %m`3/s`2

h_GOCE = 225; % km

% r and v are the initial position and velocity
[r,v] = kep2cart(I,Omega,w,M,e,a,GM); % meter for r and m/s for v
r11 = [r';v'];
f1 = drag_force(dc,h_GOCE,v');
r12 = [r';v';f1];


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


%% Aerobraking
a_Aero = (120e3+6378137+1000e3+6378137)/2; % meter
e_Aero = (1000e3-120e3)/(1000e3+120e3); % keine Einheit


h_Aero= 1000; % max altitude in km

% r and v are the initial position and velocity
[r_Aero,v_Aero] = kep2cart(I,Omega,w,M,e_Aero,a_Aero,GM); % meter for r and m/s for v
r11_Aero = [r_Aero';v_Aero'];
f1_Aero = drag_force(dc,h_Aero,v_Aero');
r12_Aero = [r_Aero';v_Aero';f1_Aero];

TC_Aero=2*pi*(sqrt((norm(r11_Aero(1:3)))^3/GM));  % aus dem Bachelor übernommen
t_1_sec_Aero=[0 10*TC_Aero];                       % aus dem Bachelor übernommen
[T1_Aero,Y1_Aero]=ode45(@odefun,t_1_sec_Aero,r12_Aero,options);

figure;
[x,y,z]=ellipsoid(0,0,0,6378137,6378137,6356752.3142);
surf(x, y, z);
axis equal;
grid on;
hold on;
title('Umlaufbahn Aerobraking')
plot3(Y1_Aero(:,1),Y1_Aero(:,2),Y1_Aero(:,3),'b');

% numerische Integration
function dydt =odefun(t,y)
GM=3.9865005e14;
dydt=[y(4);y(5);y(6);-y(1)*GM/norm(y(1:3))^3+y(7);-y(2)*GM/norm(y(1:3))^3+y(8);-y(3)*GM/norm(y(1:3))^3+y(9);0;0;0];
% Ich musste dydt durch die 0;0;0 ergänzen, da r12 9 Elemente hat und somit
% dydt auh 9 Elemente braucht. Die 9 Elemente von r12 kommen 3x Position
% 3x Geschwindigkeit und 3x Drag_force
end