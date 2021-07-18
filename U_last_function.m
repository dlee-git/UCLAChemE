function iteration = U_last_function(yBN, yUN1, yUN)
%Define Step Size 
nodes = 100; h = 1/(nodes-1);
% Define Constants
d = 0.0; Z = 0.0; D = 0.1; B = 1.5; Y = 0.05; E = 0.0; H = 0.00; theta = 0.0;

iteration = D*(2*yUN1 - 2*yUN*(1+h*theta))/h^2 + d*yBN - Z*yUN;

end