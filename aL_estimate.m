function [L atd vxi vxf vavg L_est atd_est] = aL_estimate(k, A, t_c, m, v)

%% Bounds 
%first index is L0
%second index is a_td

lb = [.8, 1.1];
ub = [1.2, 1.3];

%% Set Constants
     
    g = 9.80665; %acceration due to gravity (m/s^2)
    
   
    F_v = k*(A-g/8*t_c^2)/(1-k/m*(t_c/pi)^2); %vGRF max estimate per Morin et al. (2005)
    d_L = (A-g/8*t_c^2)/(1-k/m*(t_c/pi)^2);
    B = (pi*m*g)/(2*F_v);
    t_f = t_c*(1/B-1);
    d_y = (F_v*t_c^2)/(m*pi^2)-g/8*(t_c)^2; %vertical CoM disp. at midstance per Morin et al. (2005)
    y_f = (g/8*t_f^2);
    

%% Solve for horizontal velocties and displacement
   
syms L0

L = vpasolve(0 == ((v-(2*L0*B/t_c*(1-(L0-A)^2/L0^2)^.5))/(1-B))^2 ...
                    - (4*L0/t_c*(1-(L0-A)^2/L0^2)^.5-((v-(2*L0*B/t_c*(1-(L0-A)^2/L0^2)^.5))/(1-B)))^2 ...
                    + 2*g*(d_y+y_f) - k/m*d_L^2, L0, [0.5,1.5]);

L = double(L);
atd = asin((L-A)/L);
vxi = ((v-(2*L*B/t_c*(1-(L-A)^2/L^2)^.5))/(1-B));
vxf = (4*L/t_c*(1-(L-A)^2/L^2)^.5-((v-(2*L*B/t_c*(1-(L-A)^2/L^2)^.5))/(1-B)));
vavg = (vxi+vxf)/2;


xd = vavg*t_c/2;

if L>= lb(1) && L<= ub(1) && atd >= lb(2) && atd <= ub(2)
    L_est = L;
    atd_est = atd*180/pi;
else

%% If estimates are out of bounds, approximate L0 and a_td based on selected bounds
%Alternatively, code in use of conventional estiamte/approximation (e.g.
%atd from speed)

parameterfun = @(x, xd)(abs(xd-x(1)*cos(x(2))));
f = @(x)parameterfun(x, xd);

x0=[.9, acos(v*t_c/2)];                

options = optimoptions(@fmincon,'Algorithm', 'sqp','TolFun',1e-50, 'StepTolerance', 1e-50, 'MaxFunctionEvaluations', 200, 'Display', 'off');

x = fmincon(f,x0,[],[],[],[],lb,ub, [], options);

L_est = x(1); 
atd_est = acos(xd/L_est)*180/pi;

end

end