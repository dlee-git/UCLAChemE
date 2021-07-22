clear all; close all; clc;

% Given parameters
x = [0.40, 0.24, 0.13, 0.23];
y = 0.2;
% Guess the other y values

% Generating initial guess for temperature & Pressure
P = 0.1;
T = 100;
psat = f_Psat(T);

% Initial assumption(s): phi=1, gamma = 1
phi= [1,1,1,1];
gamma= [1,1,1,1];

% Setting loop parameters
count = 0;
toly = 0.0001;
errP = 100;
raoultsCheck = 100;
tolP = 0.005;

while abs(errP) > tolP && P < 4
    P = P + 0.01;
    T = 100;
    erry = 100;
    while abs(erry) > toly && T < 400
        if erry > 0.1
            T = T + 1;
        else
            T = T+ 0.01;
        end
        count = count+1;
        psat = f_Psat(T);
        gamma = f_gamma(T,x);
        y(2:4) = gamma(2:4).*psat(2:4).*x(2:4)./(P*phi(2:4));
        phi = f_phi(P,T,y);
        P_guess = sum(gamma.*psat.*(x)./phi);
        errP = P_guess-P;
        raoultsCheck = sum(P*y.*phi - x.*gamma.*psat);
        erry = 1 - sum(y);
    end
end

fprintf('ITERATION %.0f\n',count);
fprintf('T = %.3f\n',T);
fprintf('y1 = %.3f\n',y(1));
fprintf('y2 = %.3f\n',y(2));
fprintf('y3 = %.3f\n',y(3));
fprintf('y4 = %.3f\n',y(4));
fprintf('errorP = %.3d\n',errP);
fprintf('errorY = %.3d\n',erry);
fprintf('P = %.3f\n',P);