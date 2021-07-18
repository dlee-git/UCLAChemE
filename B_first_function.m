function iteration = B_first_function(yA1, yB1, yB2, yU1)
%Define Step Size 
nodes = 100; h = 1/(nodes-1);
% Define Constants
d = 0.0; Z = 0.0; D = 0.1; B = 1.5; Y = 0.05; E = 0.0; H = 0.00; theta = 0.0;

iteration = D*(2*yB2 - 2*yB1*(1+h) +2*h*B)/h^2 - 2*yA1*(yB1)^2 - d*yB1 + Z*yU1;

end