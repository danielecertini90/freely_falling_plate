% "﻿Analysis of transitions between fluttering, tumbling and steady descent
% of falling cards" by Andersen, Pesavento, Wang 2005

% posso avere valori super precisi di x ed y per le condizioni iniziali ed
% il resto del grafico usando web plot digitizer, non penso valga la pena
% investirci del tempo

% controllare i grafici, sono molto sensibili alle condizioni iniziali
% controllare i periodi in ogni grafico

% francesco dice di riprodurre tutti i risultati di Wang, i sei grafici (da
% (a) a (f) di Fig.3) 
% Fare anche quello con i parametri di Alsomitra che potrebbero dare il
% moto di Alsomitra che non deve ruotare, ossia deve avere theta periodico
% tipo tra -30 < theta < 30.

% supportato da Andrea, Francesco, Simone e Nicola
% clear all ; close all; clc

f=12;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');
set(0,'defaultTextInterpreter','latex');
set(0, 'defaultAxesFontSize', f)
set(0, 'defaultLegendFontSize', f)
set(0, 'defaultAxesFontName', 'Times New Roman');
set(0, 'defaultLegendFontName', 'Times New Roman');
set(0, 'DefaultLineLineWidth', 1.0);
paperUnits = 'centimeters';
paperPosition = [0 0 15 7.5];

% non dimensional period of oscillation T, Fig. 1-a
T = 8.1;

% time interval over which to solve the ODEs
tRange = [0 4*T];
% Period obtained by looking at Fig. 3-a
% tRange = [0 4.10*T];

% initial conditions for v_xp, v_yp, omega, theta, x, y
% the initial condition theta = pi/4, x~=4 and y=0  can be obtained from Fig. 3a
Y0 = [0.001; 0.001; 0.001; pi/4; 4; 0]; 

% solution of the ODE set
[tSol,ySol] = ode45(@alsomitragliding_3a,tRange,Y0);

% extract the single variables from the vector with the solutions
% x component of velocity in the reference system of the body
v_xp = ySol(:,1);

% y component of velocity in the reference system of the body
v_yp = ySol(:,2);

% omega, first derivative of theta
omega = ySol(:,3);

% theta angle defined in Fig.2
theta = ySol(:,4);

% x, horizontal coordinate in reference system linked to the lab
x_ = ySol(:,5);

% y, horizontal coordinate in reference system linked to the lab
y_ = ySol(:,6);

% transform v_xp and v_yp (velocity components in the coordinate system following the 
% rotation of the card) in v_x and v_y (velocity components in the fixed coordinate system).

v_x = v_xp.*cos(theta) - v_yp.*sin(theta);

v_y = v_xp.*sin(theta) + v_yp.*cos(theta);

figure
plot(tSol,v_xp,'k')
title('Velocity component vxp of the centre of mass, coordinate system following the rotation of the card')
xlabel('Time $t$','FontSize',f)
ylabel('$v_xp$','FontSize',f)

figure
plot(tSol,v_yp,'k')
title('Velocity component vyp of the centre of mass, coordinate system following the rotation of the card')
xlabel('Time $t$','FontSize',f)
ylabel('$v_yp$','FontSize',f)

figure
plot(tSol,v_x,'k')
title('Velocity component vx of the centre of mass, fixed coordinate system in the laboratory reference frame')
xlabel('Time $t$','FontSize',f)
ylabel('$v_x$','FontSize',f)

figure
plot(tSol,v_y,'k')
title('Velocity component vy of the centre of mass, fixed coordinate system in the laboratory reference frame')
xlabel('Time $t$','FontSize',f)
ylabel('$v_y$','FontSize',f)

figure
plot(tSol,omega,'k')
title('omega, angular velocity coordinate system following the rotation of the card')
xlabel('Time $t$','FontSize',f)
ylabel('$\omega$','FontSize',f)

figure
plot(tSol,theta,'k')
title('theta, angle in the coordinate system following the rotation of the card')
xlabel('Time $t$','FontSize',f)
ylabel('$\theta$','FontSize',f)

figure
plot(x_,y_,'k')
title('Fig. 3a Andersen, Pesavento, Wang 2005')
xlabel('$x$ horizontal displacement','FontSize',f)
ylabel('$y$ vertical displacement','FontSize',f)
axis([-10 10 -20 0])
axis equal 

% For Chandan
figure
plot(tSol,x_,'k')
title('Displacement on $x$ versus Time')
xlabel('$t$ Time','FontSize',f)
ylabel('$x$ horizontal displacement','FontSize',f)
%axis([-10 10 -20 0])
axis equal 

% For Chandan
figure
plot(tSol,y_,'k')
title('Displacement on $y$ versus Time')
xlabel('$t$ Time','FontSize',f)
ylabel('$y$ horizontal displacement','FontSize',f)
%axis([-10 10 -20 0])
axis equal 


% function dydt = alsomitragliding(t,y)
%     % Extract v_xp, v_yp, omega, theta x_ and y_ from input vector y
%     v_xp = y(1);
%     v_yp = y(2);
%     omega = y(3);
%     theta = y(4);
%     x_ = y(5);
%     y_ = y(6);
%     
%     % Define the constants Iast, C_T, C_R, A, B, mu1, mu2
%     Iast = 1.1;
%     C_T = 1.2;
%     C_R = pi;
%     A = 1.4;
%     B = 1.0;
%     mu1 = 0.2;
%     mu2 = mu1;
%     
%     % the following equations are all checked with syms Matlab, in the
%     % online course https://matlabacademy.mathworks.com/R2020a/portal.html?course=symbolic#chapter=8&lesson=1&section=1
%     % you can type the symbolic equation and get a graphical representation
%     % as in Latex.
%     
%     %%%%syms Gamma Iast pi C_T C_R A B mu1 mu2 theta omega v_xp v_yp F_vx F_vy Tau_v
%     
%     % Gamma eq. 4.7
%     Gamma = 2/pi*(-C_T*(v_xp*v_yp/sqrt(v_xp^2+v_yp^2))+C_R*omega);
%     % F_vx eq.4.8
%     F_vx = 1/pi*(A-B*((v_xp^2-v_yp^2)/(v_xp^2+v_yp^2)))*sqrt(v_xp^2+v_yp^2)*v_xp;
%     % F_vy eq.4.8
%     F_vy = 1/pi*(A-B*((v_xp^2-v_yp^2)/(v_xp^2+v_yp^2)))*sqrt(v_xp^2+v_yp^2)*v_yp;
%     % Tau eq.4.9
%     Tau_v = (mu1 + mu2*abs(omega))*omega;
%     
%     % Define dv_xpdt, dv_ypdt, domegadt, dthetadt, dx_dt and dy_dt from the ODEs
%     % ripartire da qui copiando le (5.1 - 5.2 - 5.3)
%     dv_xpdt = 1/Iast*((Iast + 1)*omega*v_yp - Gamma*v_yp - sin(theta) - F_vx);  
%     dv_ypdt = 1/(Iast + 1)*(-Iast*omega*v_xp + Gamma*v_xp - cos(theta) - F_vy); 
%     domegadt = 4/(Iast+1/2)*(-v_xp*v_yp-Tau_v); 
%     dthetadt = omega; 
%     dx_dt = v_xp.*cos(theta) - v_yp.*sin(theta);
%     dy_dt = v_xp.*sin(theta) + v_yp.*cos(theta);
%     
%      % Create output column vector dydt
%     dydt = [dv_xpdt; dv_ypdt; domegadt; dthetadt; dx_dt; dy_dt];
% end
