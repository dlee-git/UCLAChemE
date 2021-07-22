clear all; close all; clc;

% Given parameters
x = [0.40, 0.24, 0.13, 0.23];

% Generating initial guess for temperature
T = 250;
P = 0.5;
psat = f_Psat(T);

% Initial assumption(s): phi=1, gamma = 1
phi= [1,1,1,1];
gamma= [1,1,1,1];
V_water = 0.2;

% Setting loop parameters
count = 0;
errY = 100;
tolY = 0.00001;

% Iterate until water composition in vapor phase match
while abs(errY)>tolY
    % Update Temperature
    if errY > 0.1
        T = T + 1;
    else
        T = T+ 0.05;
    end
    count = count+1;
    psat = f_Psat(T);
    gamma = f_gamma(T,x);
    y = gamma.*psat.*x./(P*phi);
    phi = f_phi(P,T,y);
    P = sum(sum(gamma.*psat.*(x)./phi));
    % Calculate difference between vapor composition of water in
    % calculation vs given in the prompt
    errY = y(1) - V_water;
end


fprintf('ITERATION %.0f\nComposition of Vapor Phase, Temperature, and Pressure Calculation Real\n',count);
fprintf('T(K) = %.3f\n',T);
fprintf('y(water) = %.3f\n',y(1));
fprintf('y(ethanol) = %.3f\n',y(2));
fprintf('y(acetone) = %.3f\n',y(3));
fprintf('y(acetic acid) = %.3f\n',y(4));
fprintf('error = %.3d\n',errY);
fprintf('P(atm) = %.3f\n',P);