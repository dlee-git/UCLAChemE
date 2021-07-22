clear all; close all; clc;

% Given parameters
x = [0.40, 0.24, 0.13, 0.23];
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
    if errP > 0.1
        T = T + 1;
    else
        T = T+ 0.05;
    end
    count = count+1;
    psat = f_Psat(T);
    % Solve the composition of first bubble 
    y = gamma.*psat.*x./P;
    % Calculate the Bubble Point Pressure at given temperature
    P_guess = sum(sum(gamma.*psat.*(x)./phi));
    % Calculate difference between guess and P defined in prompt
    errP = P_guess-P;
end

% Print out results
fprintf('ITERATION %.0f\nBubble Point Calculation Ideal\n',count);
fprintf('T(K) = %.3f\n',T);
fprintf('y(water) = %.3f\n',y(1));
fprintf('y(ethanol) = %.3f\n',y(2));
fprintf('y(acetone) = %.3f\n',y(3));
fprintf('y(acetic acid) = %.3f\n',y(4));
fprintf('error = %.3d\n',errP);
fprintf('P(atm) = %.3f\n',P);
