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