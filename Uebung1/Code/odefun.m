function dydt =odefun(t,y)
GM=3.9865005e14;
dydt=[y(4);y(5);y(6);-y(1)*GM/norm(y(1:3))^3+y(7);-y(2)*GM/norm(y(1:3))^3+y(8);-y(3)*GM/norm(y(1:3))^3+y(9);0;0;0];
% Ich musste dydt durch die 0;0;0 ergänzen, da r12 9 Elemente hat und somit
% dydt auh 9 Elemente braucht. Die 9 Elemente von r12 kommen 3x Position
% 3x Geschwindigkeit und 3x Drag_force
end