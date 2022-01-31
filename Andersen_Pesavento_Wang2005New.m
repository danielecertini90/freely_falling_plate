% solution following strictly the paper
% "﻿Analysis of transitions between fluttering, tumbling and steady descent
% of falling cards" by Andersen, Pesavento, Wang 2005

clear all 
close all
clc

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

% non dimensional period of oscillation T, Fig. 1a
T = 8.1;

% time interval over which to solve the ODEs
% Maybe I should not start from 0
% tRange = [0 3*T];
tRange = [0 6*T];

% initial conditions for v_xp, v_yp, omega, theta
% the initial condition theta = pi/4 can be obtained by Fig. 3a
Y0 = [0.1; 0.1; 0.1; pi/4]; 

% solution of the ODE set
[tSol,ySol] = ode45(@alsomitragliding,tRange,Y0);

% extract the single variables from the vector with the solutions
% x component of velocity in the reference system of the body
v_xp = ySol(:,1);

% y component of velocity in the reference system of the body
v_yp = ySol(:,2);

% omega, first derivative of theta
omega = ySol(:,3);

% theta angle defined in Fig.2
theta = ySol(:,4);

% transform v_xp and v_yp (velocity components in the coordinate system following the 
% rotation of the card) in v_x and v_y (velocity components in the fixed coordinate system).

v_x = v_xp.*cos(theta) - v_yp.*sin(theta);

v_y = v_xp.*sin(theta) + v_yp.*cos(theta);


% I need to find x and y integrating v_x and v_y

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

function dydt = alsomitragliding(t,y)
    % Extract v_xp, v_yp, omega and theta from input vector y
    v_xp = y(1);
    v_yp = y(2);
    omega = y(3);
    theta = y(4);
    
    % Define the constants Iast, C_T, C_R, A, B, mu1, mu2
    Iast = 1.1;
    C_T = 1.2;
    C_R = pi;
    A = 1.4;
    B = 1.0;
    mu1 = 0.2;
    mu2 = mu1;
    
    % the following equations are all checked with syms Matlab, in the
    % online course https://matlabacademy.mathworks.com/R2020a/portal.html?course=symbolic#chapter=8&lesson=1&section=1
    % you can type the symbolic equation and get a graphical representation
    % as in Latex.
    
    %%%%syms Gamma Iast pi C_T C_R A B mu1 mu2 theta omega v_xp v_yp F_vx
    %%%%F_vy Tau_v
    
    % Gamma eq. 4.7
    Gamma = 2/pi*(-C_T*(v_xp*v_yp/sqrt(v_xp^2+v_yp^2))+C_R*omega);
    % F_vx eq.4.8
    F_vx = 1/pi*(A-B*((v_xp^2-v_yp^2)/(v_xp^2+v_yp^2)))*sqrt(v_xp^2+v_yp^2)*v_xp;
    % F_vy eq.4.8
    F_vy = 1/pi*(A-B*((v_xp^2-v_yp^2)/(v_xp^2+v_yp^2)))*sqrt(v_xp^2+v_yp^2)*v_yp;
    % Tau eq.4.9
    Tau_v = (mu1 + mu2*abs(omega))*omega;
    
    % Define dv_xpdt, dv_ypdt, domegadt and dthetadt from the ODEs
    % ripartire da qui copiando le (4.7 - 4.8 - 4.9 -5.1 - 5.2 - 5.3)
    dv_xpdt = 1/Iast*((Iast + 1)*omega*v_yp - Gamma*v_yp - sin(theta) - F_vx);  
    dv_ypdt = 1/(Iast + 1)*(-Iast*omega*v_xp + Gamma*v_xp - cos(theta) - F_vy); 
    domegadt = 4/(Iast+1/2)*(-v_xp*v_yp-Tau_v); 
    dthetadt = omega; 
    
     % Create output column vector dydt
    dydt = [dv_xpdt; dv_ypdt; domegadt; dthetadt];
end
