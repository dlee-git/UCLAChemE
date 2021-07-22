clc; clear all; close all;

 
index = 1;
for i = 1:4
    for q = 0:0.1:1
        x = zeros(4, 1);
        if i == 1   x(1) = q; else  x(1) = 1; end
        if i == 2   x(2) = q; elseif i < 2  x(2) = 0; else x(2) = 1; end
        if i == 3   x(3) = q; elseif i < 3  x(3) = 0; else x(3) = 1; end
        if i == 4   x(4) = q; elseif i < 4  x(4) = 0; else x(4) = 1; end
        soln = fsolve(@eqns, x);
        for count = 1:4
            solution(index, count) = soln(count, 1);
        end
        index = index + 1;
    end
end

function F = eqns(x)
x30 = 0.025;
x40 = 0.025;
E1 = 1.0066;
E2 = 1.0532;
U = 8;
eps1 = 12.5;
eps2 = 40;
eps3 = 1;
Da1 = 10^6;
Da2 = 10^7;
Da1a = 1.5*10^6;
Da2a = 0;
 
F(1) = 1 - x(1,1) - (Da1*exp(-E1/x(3,1))*x(1,1));
F(2) = -x(2,1) + (Da1*exp(-E1/x(3,1))*x(1,1)) - (Da2*exp(-E2/x(3,1))*x(2,1));
F(3) = x30 - x(3,1) + (Da1a*exp(-E1/x(3,1))*x(1,1)) + (Da2a*exp(-E2/x(3,1))*x(2,1)) - (U*(x(3,1) - x(4,1)));
F(4) = eps1 * eps2 * (x40 - x(4,1)) + U*eps1*eps3*(x(3,1) - x(4,1));

end