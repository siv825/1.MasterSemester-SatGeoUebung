%% Satenllitengeodäsie Übung 1 Nadine Sprügel und Toben Blei
GM=3.9865005*10^14;

%Sat 1: CHAMP
CHAMP.X0=1358287.91690107;	
CHAMP.Y0=-590375.336194285;
CHAMP.Z0=-6632649.79730091;
CHAMP.VX0=-6036.56551965113;
CHAMP.VY0=4447.8030057021;
CHAMP.VZ0=-1632.11801174039;
yC0=[CHAMP.X0;CHAMP.Y0;CHAMP.Z0;CHAMP.VX0;CHAMP.VY0;CHAMP.VZ0];

%Sat 2: Molnija
Molnija.X0=2821874.45574537;	
Molnija.Y0=-1315861.66906604;	
Molnija.Z0=-6110772.2247106;
Molnija.VX0=4249.95368432902;	
Molnija.VY0=9114.05508788868;
Molnija.VZ0=-1.64595869549316e-12;
yM0=[Molnija.X0; Molnija.Y0;Molnija.Z0;Molnija.VX0;Molnija.VY0;Molnija.VZ0];

% Umlaufzeit CHAMP
eC=0;
TC=2*pi*(sqrt((norm(yC0(1:3)))^3/GM));
NC=2*pi/TC;

%Umlaufzeit Molnija
aM=26378000;
TM=2*pi*(sqrt((aM^3)/GM));

%% numerische Integration
options = odeset('RelTol',1e-10,'AbsTol',1e-15);
%CHAMP
tspanC=[0 5*TC];
[tC,yC]=ode45(@odefun,tspanC,yC0,options);


%Molnija
tspanM=[0 5*TM];
[tM,yM1]=ode45(@odefun,tspanM,yM0);


%% Plot

figure;
[x,y,z]=ellipsoid(0,0,0,6378137,6378137,6356752.3142);
surf(x, y, z);
axis equal;
grid on;
hold on;
title('Umlaufbahn CHAMP')
plot3(yC(:,1),yC(:,2),yC(:,3),'LineWidth',2);



figure;
[x,y,z]=ellipsoid(0,0,0,6378137,6378137,6356752.3142);
surf(x, y, z);
axis equal;
grid on;
hold on;
title('Umlaufbahn Molnija')
plot3(yM(:,1),yM(:,2),yM(:,3),'LineWidth',2);

figure;
[x,y,z]=ellipsoid(0,0,0,6378137,6378137,6356752.3142);
surf(x, y, z);
axis equal;
grid on;
hold on;
title('Umlaufbahn CHAMP und Molnija')
plot3(yC(:,1),yC(:,2),yC(:,3),'b','LineWidth',2);
hold on;
plot3(yM(:,1),yM(:,2),yM(:,3),'r','LineWidth',2);

%% noch berechnen: Energie des Satelliten, Flächengeschwindigkeit C, Winkel zwischen der momentalen Bahnebene und der Anfangsbahnebene

% Energie des Satelliten
for i=1:length(yC)
    ECt(i)=0.5*(norm(yC(i,4:6))^2)-GM/(norm(yC(i,1:3)));
end
    EC=ECt';
    figure;
    plot(EC);
    title('Energie des Satellitens CHAMP')
   
    
for i=1:length(yM)
    EMt(i)=0.5*(norm(yM(i,4:6))^2)-GM/(norm(yM(i,1:3)));
end
    EM=EMt';
    figure;
    plot(EM);
    title('Energie des Satellitens Molnija')
    
% Flächengeschwindigkeit C
kC=cross(yC(:,1:3),yC(:,4:6));
CC=0.5*vecnorm(kC');
figure;
plot(CC);
title('Flächengeschwindigkeit CHAMP')

kM=cross(yM(:,1:3),yM(:,4:6));
CM=0.5*vecnorm(kM');
figure;
plot(CM);
title('Flächengeschwindigkeit Molnija')


% Winkel zwischen der momentanen Bahnebene und der Anfangsbahnebene

LCHAMP = cross(yC0(1:3),yC0(4:6));
LMOLNIJA = cross(yM0(1:3),yM0(4:6));


% Kleinwinkelnährung verwendet sin(x)~x
for i=1:length(kC)
    phiC(i) = norm(cross(kC(i,1:3),LCHAMP))/(norm(LCHAMP)*norm(kC(i,1:3)));
end
figure
plot(phiC)
title('Winkel zwischen der momentanen Bahnebene und der Anfangsbahnebene des Satellitens CHAMP')


for i=1:length(kM)
    phiM(i) = norm(cross(kM(i,1:3),LMOLNIJA))/(norm(LMOLNIJA)*norm(kM(i,1:3)));
end
figure
plot(phiM)
title('Winkel zwischen der momentanen Bahnebene und der Anfangsbahnebene des Satellitens Molnija')




% numerische Integration
function dydt =odefun(t,y)
GM=3.9865005e14;
dydt=[y(4);y(5);y(6);-y(1)*GM/norm(y(1:3))^3;-y(2)*GM/norm(y(1:3))^3;-y(3)*GM/norm(y(1:3))^3];
end








