clear all; close all; clc;
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
tic();

% non dimensional period of oscillation T, Fig. 1-b
T = 10;

% time interval over which to solve the ODEs
tRange = [0 T];

% initial conditions for u, v, omega, theta, (x, y), theta is 0
% the initial condition theta = -pi/4, x=0 and y=0  can be obtained from Fig.3 
Y0 = [0.1; 0.1; 0.1; pi/4; 0.1; 0.1]; 

% solution of the ODE set
[tSol,ySol] = ode45(@Huang_freelyfallingplate,tRange,Y0);

% extract the single variables from the vector with the solutions
% x component of velocity in the reference system of the body
u = ySol(:,1);

% y component of velocity in the reference system of the body
v = ySol(:,2);

% omega, first derivative of theta
omega = ySol(:,3);

% theta angle defined in Fig.2
theta = ySol(:,4);

% x, horizontal coordinate in reference system linked to the lab
x_ = ySol(:,5);
 
% y, horizontal coordinate in reference system linked to the lab
y_ = ySol(:,6);

figure
plot(tSol,u,'k')
title('Velocity component u of the centre of mass, coordinate system following the rotation of the card')
xlabel('Time $t$','FontSize',f)
ylabel('$u$','FontSize',f)

figure
plot(tSol,v,'k')
title('Velocity component v of the centre of mass, coordinate system following the rotation of the card')
xlabel('Time $t$','FontSize',f)
ylabel('$v$','FontSize',f)

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
title('Fig. 3 Huang 2013')
xlabel('$x$ horizontal displacement','FontSize',f)
ylabel('$y$ vertical displacement','FontSize',f)
% axis([-3 57 -40 0])
axis equal 
toc();
function dydt = Huang_freelyfallingplate(t,y)
    % Extract u, v, omega, theta, (x_ and y_) from input vector y
    u = y(1);
    v = y(2);
    omega = y(3);
    theta = y(4);
    x_ = y(5);
    y_ = y(6);
    
    % Value from table 1 for a,b,c,d e/a in [cm] so /100
    a = 2.5;
%     a = 3.0;

    b = 0.3;

    c = 1;
%     c = 1.2;

    d = 0.1;

    % Fig.3 top, first from the top
    % e/a = 0, 
    e = 0;
    % Fig.3, second from the top 
%     e = 0.0086*a;
    % Fig.3, third from the top
%     e = 0.014*a;
    % Fig.3, fourth from the top, last one
%     e = 0.017*a;
    
    % Value from table 2 for C_T, C_R, C_0, C_pi, C_tau 
    C_T = 4.5;
    C_R = 1.8;
    C_0 = 0.2;
    C_pi = 0.5;
    C_tau = 1.9;
    
    % Values for rho in g/cm^3 so *1000, near eq 1 that describes m
    rho_f = 1.0;
    rho_p = 1.2;
    rho_a = 2.7;
    
    % it uses cm, so g is cm/s^2
    g = 9.80665*100;
    
    % Define the constants
    m = rho_p*a*b + (rho_a - rho_p)*c*d;
    m11 = 3*pi*rho_f*b^2/8;
    m22 = 3*pi*rho_f*a^2/8;
    mp = m - rho_f*a*b;
    
    % if e is 0, h is 0
    h = m*e/((rho_a-rho_p)*c*d);
    
    I = rho_p*a*b*(((a^2+b^2)/12)+e^2) + (rho_a-rho_p)*c*d*(((c^2+d^2)/12)+(h-e)^2);    
    Ia = (5*pi/256)*rho_f*(a^2-b^2)^2 + (3*pi/8)*rho_f*a^2*e^2;
     
    up = u;
    vp = v - e*omega;
   
    % Gamma eq.9
    Gamma = -C_T*a*up*vp/sqrt(up^2 + vp^2) + 1/2 * C_R*a^2*omega;
    
    % F_x eq.10
    F_x = (rho_f*a/2)*(C_0*up^2/sqrt(up^2+vp^2) + C_pi*vp^2/sqrt(up^2+vp^2))*up;
    
    % F_y eq.11
    F_y = (rho_f*a/2)*(C_0*up^2/sqrt(up^2+vp^2) + C_pi*vp^2/sqrt(up^2+vp^2))*vp;
    
    % Tau eq.12
    Tau = C_tau*rho_f*a^4*omega*abs(omega)/64;

    % Define dudt, dvdt, domegadt, dthetadt, dxdt and dydt from the ODEs
    % 5-6-7
    
    dudt = 1/(m + m11)*((m + m22)*omega*v - mp*g*sin(theta) - rho_f*Gamma*v - F_x);  
    
    dvdt = 1/(m + m22)*(-(m + m11)*omega*u - mp*g*cos(theta) + rho_f*Gamma*u - F_y); 
    
    domegadt = 1/(I + Ia)*((m11 - m22)*u*v - rho_f*a*b*g*e*cos(theta) - Tau);
    
    dthetadt = omega; 
    
    dx_dt = u*cos(theta) - v*sin(theta);
     
    dy_dt = u*sin(theta) + v*cos(theta);
     
    % Create output column vector dydt
    dydt = [dudt; dvdt; domegadt; dthetadt; dx_dt; dy_dt];
end