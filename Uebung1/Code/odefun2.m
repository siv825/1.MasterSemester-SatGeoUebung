
function dydt2 =odefun2(t,y,dc)
GM=3.9865005e14;
R = 6371000; % Erde Radius
h = norm(y(1:3)) - R; % Hoehe
f1_Aero = drag_force(dc,h,[y(4);y(5);y(6)]);
dydt2=[y(4);
      y(5);
      y(6);
      -y(1)*GM/norm(y(1:3))^3 + f1_Aero(1);
      -y(2)*GM/norm(y(1:3))^3 + f1_Aero(2);
      -y(3)*GM/norm(y(1:3))^3 + f1_Aero(3)];
% Ich musste dydt durch die 0;0;0 ergänzen, da r12 9 Elemente hat und somit
% dydt auh 9 Elemente braucht. Die 9 Elemente von r12 kommen 3x Position
% 3x Geschwindigkeit und 3x Drag_force
end