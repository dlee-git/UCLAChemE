clc; clear all; close all;
 

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


syms x1 x2 x3 x4 x30 x40
F = [1 - x1 - (Da1*exp(-E1/x3)*x1)
    -x2 + (Da1*exp(-E1/x3)*x1) - (Da2*exp(-E2/x3)*x2)
    x30 - x3 + (Da1a*exp(-E1/x3)*x1) + (Da2a*exp(-E2/x3)*x2) - (U*(x3 - x4))
    eps1 * eps2 * (x40 - x4) + U*eps1*eps3*(x3 - x4)];
v = [x1, x2, x3, x4, x30, x40];
 
jac = jacobian(F, v);
 
linear = double(subs(jac, v, [0.7877 0.0908 0.0665 0.0319 0.025 0.025]));
