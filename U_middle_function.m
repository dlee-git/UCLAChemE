function iteration = U_middle_function(yB2, yU1, yU2, yU3)
%Define Step Size 
nodes = 100; h = 1/(nodes-1);
% Define Constants
d = 0.0; Z = 0.0; D = 0.1; B = 1.5; Y = 0.05; E = 0.0; H = 0.00; theta = 0.0;

iteration = D*(yU1 - 2*yU2 +yU3)/h^2 + d*yB2 - Z*yU2;

end