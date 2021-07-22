clear all; close all; clc;

% Given parameters
T = 273;
P = 1;
x = [40, 24, 13, 23];

% Initial assumption(s): phi=1
phi= [1;1];

% Generating initial guess for pressure
psat = f_Psat(T);
gamma = f_gamma(T,x);
% P = sum(gamma.*psat.*(x)./phi);

% Setting loop parameters
c = 0;
errP = 100;
tolP = 10^-4;

% Loop shown in Figure (14.1)
while abs(errP)>tolP
    fprintf('------ Press any key to continue ------\n')
    pause;
    c = c+1;
    T = T + 1;
    gamma = f_gamma(T,x);
    psat = f_Psat(T);
    y = gamma.*psat.*x/P;
    P_guess = sum(gamma.*psat.*(x)./phi);
    phi = f_phi(P,T,y);
    errP = P_guess-P;
    fprintf('ITERATION %.0f\n',c);
    fprintf('P = %.3f\n',P);
    fprintf('y1 = %.3f\n',y(1));
    fprintf('y2 = %.3f\n',y(2));
    fprintf('error = %.3d\n',errP);
end


fprintf('T = %.3f\n',T);