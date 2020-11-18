%% Satgeo Ü1
% Nadine Sprügel 3317570
% Ziqing Yu 3218051

load('rhocoe.mat')
%% GOCE
a = 6378137+225e3; % meter
I = deg2rad(96.6); % radiant
e = 0; % keine Einheit
Omega = deg2rad(335);% radiant
w = deg2rad(273); % radiant
M = deg2rad(5); % radiant
GM = 3.9865005*10^14; %m`3/s`2

h_GOCE = 225 * 1000; % km

% r and v are the initial position and velocity
[r,v] = kep2cart(I,Omega,w,M,e,a,GM); % meter for r and m/s for v
r11 = [r';v'];
f1 = drag_force(dc,h_GOCE,v');



options = odeset('RelTol',1e-15,'AbsTol',1e-15);

TC=2*pi*sqrt(a^3/GM);  % aus dem Bachelor übernommen einheit:Sekunde
t_1_sec=[0 8*TC];                       % aus dem Bachelor übernommen
[T1,Y1]=ode45(@(t,y)odefun(t,y,dc),t_1_sec,r11,options);


% weil die Erde dreht sich
Y1e = zeros(length(T1),3);
theta_gr=2*pi/(24*3600)*T1;
for i=1:length(T1)
    Y1e(i,:)=R3(theta_gr(i))*Y1(i,1:3)';
end

figure;
[x,y,z]=ellipsoid(0,0,0,6378137,6378137,6356752.3142);
surf(x, y, z);
axis equal;
grid on;
hold on;
title('Umlaufbahn GOCE')
plot3(Y1e(:,1),Y1e(:,2),Y1e(:,3),'b','LineWidth',2);



%% Aerobraking
a_Aero = (120e3+6378137+1000e3+6378137)/2; % meter
% e_Aero = (1000e3-120e3)/(1000e3+120e3); % keine Einheit
e_Aero = (a_Aero - 120e3 - 6378137)/a_Aero; % keine Einheit

% r and v are the initial position and velocity
[r_Aero,v_Aero] = kep2cart(I,Omega,w,M,e_Aero,a_Aero,GM); % meter for r and m/s for v
r11_Aero = [r_Aero';v_Aero'];

% 
TC_Aero=2*pi*sqrt(a_Aero^3/GM);  % aus dem Bachelor übernommen
t_1_sec_Aero=[0 TC_Aero*10];                       % aus dem Bachelor übernommen
[T1_Aero,Y1_Aero]=ode45(@(t,y)odefun2(t,y,dc),t_1_sec_Aero,r11_Aero,options);

Y2e = zeros(length(T1_Aero),3);
theta_gr=2*pi/(24*3600)*T1_Aero;
for i=1:length(T1)
    Y2e(i,:)=R3(theta_gr(i))*Y1_Aero(i,1:3)';
end

figure;
[x,y,z]=ellipsoid(0,0,0,6378137,6378137,6356752.3142);
surf(x, y, z);
axis equal;
grid on;
hold on;
title('Umlaufbahn Aerobraking')
plot3(Y2e(:,1),Y2e(:,2),Y2e(:,3),'b','LineWidth',1.5);

