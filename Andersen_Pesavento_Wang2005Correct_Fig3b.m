% "Analysis of transitions between fluttering, tumbling and steady descent
% of falling cards" by Andersen, Pesavento, Wang 2005

clear all
close all
clc

f=16;
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

% non dimensional period of oscillation T, Fig. 1-b
T = 17.3;

% time interval over which to solve the ODEs
% Maybe I should not start from 0? It looks starting from 0 it is fine
% tRange = [0 3*T];
tRange = [0 4*T];

% initial conditions for v_xp, v_yp, omega, theta, x, y
% the initial condition theta = pi/4, x~=4 and y=0  can be obtained from Fig. 3a
Y0 = [0.001; 0.001; 0.001; 4*pi/9; 80; -300];

% solution of the ODE set
[tSol,ySol] = ode45(@alsomitragliding_3b,tRange,Y0);

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
%title('Velocity component vx of the centre of mass, fixed coordinate system in the laboratory reference frame')
xlabel('Time $t$','FontSize',f)
ylabel('$v_x$','FontSize',f)

figure
plot(tSol,v_y,'k')
%title('Velocity component vy of the centre of mass, fixed coordinate system in the laboratory reference frame')
xlabel('Time $t$','FontSize',f)
ylabel('$v_y$','FontSize',f)

figure
plot(tSol,omega,'k')
%title('omega, angular velocity coordinate system following the rotation of the card')
xlabel('Time $t$','FontSize',f)
ylabel('$\omega$','FontSize',f)

figure
plot(tSol,theta,'k')
%title('theta, angle in the coordinate system following the rotation of the card')
xlabel('Time $t$','FontSize',f)
ylabel('$\theta$','FontSize',f)

figure
plot(v_x,v_y,'k')
%title('Behavior of a Falling Paper')
xlabel('$v_x$','FontSize',f)
ylabel('$v_y$','FontSize',f)

figure
plot(x_,y_,'k')
%title('Fig. 3b Andersen, Pesavento, Wang 2005')
xlabel('$x$ horizontal displacement','FontSize',f)
ylabel('$y$ vertical displacement','FontSize',f)
% axis([-10 10 0 -20])
axis equal