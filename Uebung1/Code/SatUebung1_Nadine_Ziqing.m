%% Satellitengeodäsie MSc - Übung 1: Atmosphärenreibung

% Nadine Sprügel 3317570
% Ziqing Yu 3218051

clearvars
close all
clc
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

options = odeset('RelTol',1e-15,'AbsTol',1e-15);

TC=2*pi*sqrt(a^3/GM);  % Sekunde
t_1_sec=[0 8*TC];                       
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

% Kepler Elemente 
for i=1:length(Y1)
    [I_Gk(i),Omega_Gk(i),w_Gk(i),M_Gk(i),e_Gk(i),a_Gk(i)] = cart2kep(Y1(i,1:3)',Y1(i,4:6)',GM);
    f1_G = drag_force(dc,h_GOCE,[Y1(i,4);Y1(i,5);Y1(i,6)]);
    f1_GOCE(i,1:3)=f1_G'; %f1_GOCE=[f1, f2, f3]    
end

% Gauss LPE
for i=1:length(Y1)
    E0=M_Gk(i);
    Ek=M_Gk(i)+e_Gk(i).*sin(E0);
    while abs(E0-Ek)>1e-12
        E0=Ek;
        Ek=M_Gk(i)+e_Gk(i).*sin(E0);
    end
    E_Gk(i)=Ek;
    n(i)=sqrt(GM/a_Gk(i)^3);
    nu(i)=atan((sqrt(1-e_Gk(i)^2)*sin(E_Gk(i)))/(cos(E_Gk(i))-e_Gk(i)));
end

aPunkt=2./n'.*f1_GOCE(:,1);
ePunkt=1./(n'.*a_Gk').*2.*cos(nu').*f1_GOCE(:,1);
u=w_Gk'+nu';
%&OmegaPunkt=1/(n'.*a'.*sin(I_Gk')).*sin(u').*f1_GOCE(:,2); % Laut Folien von Omid müssen wir das nicht berechnen, brauche aber OmegaPunkt für die folgenden Formel 
OmegaPunkt=0;
omegaPunkt_MPunkt=n'-cos(I_Gk').*OmegaPunkt;


tdiff=diff(T1(1:2));

figure;
plot(T1,aPunkt)
hold on;
plot(T1(1:end-1),diff(a_Gk)./tdiff)
legend('aus Gauss LPE','aus Keplerelemente')
title('Vergleich von $\dot{a}$ (GOCE)','FontSize',15,'Interpreter','latex')

figure;
plot(ePunkt)
hold on;
plot(T1(1:end-1),diff(e_Gk)./tdiff)
legend('aus Gauss LPE','aus Keplerelemente')
title('Vergleich von $\dot{e}$ (GOCE)','FontSize',15,'Interpreter','latex')


figure;
plot(omegaPunkt_MPunkt)
title('GOCE LPE Gauss: omegaPunkt+MPunkt')
hold on;
omegaPunkt_MPunkt_Kepler=diff(w_Gk)./tdiff+diff(M_Gk)./tdiff; % woher kommen die krassen Ausschläge?
plot(T1(1:end-1),omegaPunkt_MPunkt_Kepler)
legend('aus Gauss LPE','aus Keplerelemente')
title('Vergleich von $\dot{w}+\dot{M}$ (GOCE)','FontSize',15,'Interpreter','latex')


% %% Aerobraking
% a_Aero = (120e3+6378137+1000e3+6378137)/2; % meter
% a_Aero_max = 1000e3+6378137; % meter
% a_Aero_min = (120e3+6378137); % meter
% e_Aero = (a_Aero_max-a_Aero_min)/(a_Aero_max+a_Aero_min); % keine Einheit
% 
% % r and v are the initial position and velocity
% [r_Aero,v_Aero] = kep2cart(I,Omega,w,M,e_Aero,a_Aero,GM); % meter for r and m/s for v
% r11_Aero = [r_Aero';v_Aero'];
% 
% 
% TC_Aero=2*pi*sqrt(a_Aero^3/GM);  % Sekunde
% t_1_sec_Aero=[0 TC_Aero*10];                    
% [T1_Aero,Y1_Aero]=ode45(@(t,y)odefun2(t,y,dc),t_1_sec_Aero,r11_Aero,options);
% 
% 
% % weil die Erde dreht sich
% Y2e = zeros(length(T1_Aero),3);
% theta_gr=2*pi/(24*3600)*T1_Aero;
% for i=1:length(T1_Aero)
%     Y2e(i,:)=R3(theta_gr(i))*Y1_Aero(i,1:3)';
% end
% 
% figure;
% [x,y,z]=ellipsoid(0,0,0,6378137,6378137,6356752.3142);
% surf(x, y, z);
% axis equal;
% grid on;
% hold on;
% title('Umlaufbahn Aerobraking')
% plot3(Y2e(:,1),Y2e(:,2),Y2e(:,3),'b','LineWidth',1.5);
% 
% % Kepler Elemente 
% for i=1:length(Y1_Aero)
%     [I_Ak(i),Omega_Ak(i),w_Ak(i),M_Ak(i),e_Ak(i),a_Ak(i)] = cart2kep(Y1_Aero(i,1:3)',Y1_Aero(i,4:6)',GM);
%     h = norm(Y1_Aero(i,1:3)) - 6371000;
%     f1_A = drag_force(dc,h,[Y1_Aero(i,4);Y1_Aero(i,5);Y1_Aero(i,6)]);
%     f1_Aero(i,1:3)=f1_A'; %f1_Aero=[f1, f2, f3]  
% end
% 
% % Gauss LPE
% for i=1:length(Y1_Aero)
%     E0=M_Ak(i);
%     Ek=M_Ak(i)+e_Ak(i).*sin(E0);
%     while abs(E0-Ek)>1e-12
%         E0=Ek;
%         Ek=M_Ak(i)+e_Ak(i).*sin(E0);
%     end
%     E_Ak(i)=Ek;
%     n_Aero(i)=sqrt(GM/a_Ak(i)^3);
%     nu_Aero(i)=atan((sqrt(1-e_Ak(i)^2)*sin(E_Ak(i)))/(cos(E_Ak(i))-e_Ak(i)));
% end
% 
% aPunkt_Aero=2./n_Aero'.*f1_Aero(:,1);
% ePunkt_Aero=1./(n_Aero'.*a_Ak').*2.*cos(nu_Aero').*f1_Aero(:,1);
% u_Aero=w_Ak'+nu_Aero';
% %&OmegaPunkt=1/(n'.*a'.*sin(I_Gk')).*sin(u').*f1_GOCE(:,2); % Laut Folien von Omid müssen wir das nicht berechnen, brauche aber OmegaPunkt für die folgenden Formel 
% OmegaPunkt=0;
% omegaPunkt_MPunkt_Aero=n_Aero'-cos(I_Ak').*OmegaPunkt;
% 
% 
% tdiff=diff(T1_Aero(1:2));
% 
% figure;
% plot(T1_Aero,aPunkt_Aero)
% hold on;
% plot(T1_Aero(1:end-1),diff(a_Ak)./tdiff)
% legend('aus Gauss LPE','aus Keplerelemente')
% title('Vergleich von $\dot{a}$ (Aerobrake)','FontSize',15,'Interpreter','latex')
% 
% 
% figure;
% plot(ePunkt_Aero)
% hold on;
% plot(T1_Aero(1:end-1),diff(e_Ak)./tdiff)
% legend('aus Gauss LPE','aus Keplerelemente')
% title('Vergleich von $\dot{e}$ (Aerobrake)','FontSize',15,'Interpreter','latex')
% 
% 
% figure;
% plot(omegaPunkt_MPunkt_Aero)
% hold on;
% omegaPunkt_MPunkt_Kepler_A=diff(w_Ak)./tdiff+diff(M_Ak)./tdiff; % woher kommen die krassen Ausschläge?
% plot(T1_Aero(1:end-1),omegaPunkt_MPunkt_Kepler_A)
% legend('aus Gauss LPE','aus Keplerelemente')
% title('Vergleich von $\dot{w}+\dot{M}$ (Aerobrake)','FontSize',15,'Interpreter','latex')





