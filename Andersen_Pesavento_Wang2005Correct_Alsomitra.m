% "﻿Analysis of transitions between fluttering, tumbling and steady descent
% of falling cards" by Andersen, Pesavento, Wang 2005

% posso avere valori super precisi di x ed y per le condizioni iniziali ed
% il restodel grafico usando web plot digitizer, non penso valga la pena
% investirci del tempo

% controllare i grafici, sono molto sensibili alle condizioni iniziali
% controllare i periodi in ogni grafico

% francesco dice di riprodurre tutti i risultati di Wang, i sei grafici (da
% (a) a (f) di Fig.3) 
% Fare anche quello con i parametri di Alsomitra che potrebbero dare il
% moto di Alsomitra che non deve ruotare, ossia deve avere theta periodico
% tipo tra -30 < theta < 30.

% supportato da Andrea, Francesco, Simone e Nicola
clear all 
close all
clc

%f=22;
f = 12;
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
% Periodo della dinamica di lungo periodo per la fugoide, report Ignazio
% Overleaf, T = 2*pi/omega_n*sqrt(1-chi^2)
T = 0.7;

% time interval over which to solve the ODEs
tRange = [0 100*T];
% Period obtained by looking at Fig. 3-a
% tRange = [0 4.10*T];

% initial conditions for v_xp, v_yp, omega, theta, x, y
% the initial condition theta = pi/4, x~=4 and y=0  can be obtained from Fig. 3a
Y0 = [0.001; 0.001; 0.001; pi/4; 0; 0]; 

% solution of the ODE set
[tSol,ySol] = ode45(@alsomitragliding_Alsomitra,tRange,Y0);

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
% title('Fig. 3a Andersen, Pesavento, Wang 2005')
xlabel('$x$ horizontal displacement','FontSize',f)
ylabel('$y$ vertical displacement','FontSize',f)
%axis([-10 10 -20 0])
axis equal