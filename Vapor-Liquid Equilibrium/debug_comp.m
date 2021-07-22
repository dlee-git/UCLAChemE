clear all; close all; clc;

% Given parameters
x = [0.40, 0.24, 0.13, 0.23];
y(1) = 0.2;
% Guess the other y values

% Generating initial guess for temperature & Pressure
P = 0.1;
T = sum(x.*f_Tsat(P));
psat = f_Psat(T);

% Initial assumption(s): phi=1, gamma = 1
phi= [1,1,1,1];
gamma= [1,1,1,1];

% Setting loop parameters
count = 0;
errP = 100;
tolP = 0.001;

while abs(errP)>tolP
    count = count+1;
    T = sum(x.*f_Tsat(P));
    psat = f_Psat(T);
    gamma = f_gamma(T,x);
    y(2:4) = gamma(2:4).*psat(2:4).*x(2:4)./(P*phi(2:4));
    P_guess = sum(sum(gamma.*psat.*(x)./phi));
    phi = f_phi(P,T,y);
    errP = sum(P*y.*phi - x.*gamma.*psat);
    if errP > 0.1
        P = P + 0.01;
    else
        P = P + 0.005;
    end
end

fprintf('ITERATION %.0f\n',count);
fprintf('T = %.3f\n',T);
fprintf('y1 = %.3f\n',y(1));
fprintf('y2 = %.3f\n',y(2));
fprintf('y3 = %.3f\n',y(3));
fprintf('y4 = %.3f\n',y(4));
fprintf('error = %.3d\n',errP);
fprintf('P = %.3f\n',P);