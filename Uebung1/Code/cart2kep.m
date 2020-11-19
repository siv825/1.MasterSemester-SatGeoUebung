function [i,Omega,w,M,e,a] = cart2kep(Position,Geschwindigkeit,GM)
% Berechnung der Keplerelemente aus den kartesischen Koordinaten


L=cross(Position,Geschwindigkeit);
Omega=atan2(L(1), -L(2));
i=atan2(sqrt(L(1)^2+L(2)^2),L(3));

a=GM*norm(Position)/(2*GM-norm(Position)*norm(Geschwindigkeit)^2);

if 1-(norm(L)^2/(GM*a))<0
    e=0;
else
    e=sqrt(1- norm(L)^2/(GM*a));
end

rPunkt=dot(Position,Geschwindigkeit)/norm(Position);
Ecos=(a-norm(Position))/(a*e);
Esin=norm(Position)*rPunkt/(e*sqrt(GM*a));

nu=atan2(sqrt(1-e^2)*Esin,Ecos-e);

E=atan2(Esin,Ecos);
M=E-e*sin(E);


p=R1(i)*R3(Omega)*Position;
w=atan2(p(2),p(1))-nu;
end

