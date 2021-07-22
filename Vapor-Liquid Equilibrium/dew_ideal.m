clear all; close all; clc;

% Given parameters
y = [0.40, 0.24, 0.13, 0.23];
P = 1;

% Initial assumption(s): phi=1, gamma = 1
phi= [1,1,1,1];
gamma= [1,1,1,1];

% Generating initial guess for Temperature & Saturated Pressure
T = 0;
psat = f_Psat(T);

% Setting loop parameters
count = 0;
errP = 100;
tolP = 0.001;

while abs(errP)>tolP
    % Update Temperature
    if errP > 0.1
        T = T + 1;
    else
        T = T+ 0.05;
    end
    count = count+1;
    psat = f_Psat(T);
    % Solve the composition of first liquid droplet 
    x = P*phi.*y./(gamma.*psat);
    % Calculate the Dew Point Pressure at given temperature
    P_guess = 1/(sum(phi.*y./(gamma.*psat)));
    % Calculate difference between guess and P defined in prompt
    errP = P_guess-P;
end

fprintf('ITERATION %.0f\nDew Point Calculation Ideal\n',count);
fprintf('T(K) = %.3f\n',T);
fprintf('x(water) = %.3f\n',x(1));
fprintf('x(ethanol) = %.3f\n',x(2));
fprintf('x(acetone) = %.3f\n',x(3));
fprintf('x(acetic acid) = %.3f\n',x(4));
fprintf('error = %.3d\n',errP);
fprintf('P(atm) = %.3f\n',P);