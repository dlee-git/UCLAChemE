function iteration = A_last_function(yAN1, yAN, yBN)
%Define Step Size 
nodes = 100; h = 1/(nodes-1);
% Define Constants
d = 0.0; Z = 0.0; D = 0.1; B = 1.5; Y = 0.05; E = 0.0; H = 0.00; theta = 0.0;

iteration = D*(2*yAN1 - 2*yAN*(1+h*E))/h^2 - yAN*(yBN)^2 - Y*yAN;

end