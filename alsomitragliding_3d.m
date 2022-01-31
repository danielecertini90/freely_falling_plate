function dydt = alsomitragliding_3d(t,y)
    % Extract v_xp, v_yp, omega, theta, x_ and y_ from input vector y
    v_xp = y(1);
    v_yp = y(2);
    omega = y(3);
    theta = y(4);
    x_ = y(5);
    y_ = y(6);
    
    % Define the constants Iast, C_T, C_R, A, B, mu1, mu2
    Iast = 1.6;
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
    
    %%%%syms Gamma Iast pi C_T C_R A B mu1 mu2 theta omega v_xp v_yp F_vx F_vy Tau_v
    
    % Gamma eq. 4.7
    Gamma = 2/pi*(-C_T*(v_xp*v_yp/sqrt(v_xp^2+v_yp^2))+C_R*omega);
    % F_vx eq.4.8
    F_vx = 1/pi*(A-B*((v_xp^2-v_yp^2)/(v_xp^2+v_yp^2)))*sqrt(v_xp^2+v_yp^2)*v_xp;
    % F_vy eq.4.8
    F_vy = 1/pi*(A-B*((v_xp^2-v_yp^2)/(v_xp^2+v_yp^2)))*sqrt(v_xp^2+v_yp^2)*v_yp;
    % Tau eq.4.9
    Tau_v = (mu1 + mu2*abs(omega))*omega;
    
    % Define dv_xpdt, dv_ypdt, domegadt, dthetadt, dx_dt and dy_dt from the ODEs
    % ripartire da qui copiando le (5.1 - 5.2 - 5.3)
    dv_xpdt = 1/Iast*((Iast + 1)*omega*v_yp - Gamma*v_yp - sin(theta) - F_vx);  
    dv_ypdt = 1/(Iast + 1)*(-Iast*omega*v_xp + Gamma*v_xp - cos(theta) - F_vy); 
    domegadt = 4/(Iast+1/2)*(-v_xp*v_yp-Tau_v); 
    dthetadt = omega; 
    dx_dt = v_xp.*cos(theta) - v_yp.*sin(theta);
    dy_dt = v_xp.*sin(theta) + v_yp.*cos(theta);
    
     % Create output column vector dydt
    dydt = [dv_xpdt; dv_ypdt; domegadt; dthetadt; dx_dt; dy_dt];
end



