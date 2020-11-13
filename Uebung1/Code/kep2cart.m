function [PositionT,GeschwindigkeitT] = kep2cart(i,Omega,w,Mt,e,a)
% Berechnung die kartesische Koordinaten aus Keplerelemente
% alle Winkel in rad
GM=3.9865005*10^14;
E0=Mt;
Ek=Mt+e.*sin(E0);
while abs(E0-Ek)>1e-12
    E0=Ek;
    Ek=Mt+e.*sin(E0);
end
E=Ek;
n=sqrt(GM/a^3);

for z=1:length(E)
rf=[a*(cos(E(z))-e);a*sqrt(1-e^2)*sin(E(z));0*E(z)];
vf=n*a/(1-e*cos(E(z)))*[-sin(E(z)); sqrt(1-e^2)*cos(E(z));0];
Position(z,:)=R3(-Omega)*R1(-i)*R3(-w)*rf;
Geschwindigkeit(z,:)=R3(-Omega)*R1(-i)*R3(-w)*vf;
end

PositionT=Position;
GeschwindigkeitT=Geschwindigkeit;
end

