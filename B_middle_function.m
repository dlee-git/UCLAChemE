function iteration = B_middle_function(yA2, yB1, yB2, yB3, yU2)
%Define Step Size 
nodes = 100; h = 1/(nodes-1);
% Define Constants
d = 0.0; Z = 0.0; D = 0.1; B = 1.5; Y = 0.05; E = 0.0; H = 0.00; theta = 0.0;

iteration = D*(yB1 - 2*yB2 + yB3)/h^2 - 2*yA2*(yB2)^2 - d*yB2 + Z*yU2;

end