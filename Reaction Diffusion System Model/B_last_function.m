function iteration = B_last_function(yAN, yBN1, yBN, yUN)
%Define Step Size 
nodes = 100; h = 1/(nodes-1);
% Define Constants
d = 0.0; Z = 0.0; D = 0.1; B = 1.5; Y = 0.05; E = 0.0; H = 0.00; theta = 0.0;

iteration = D*(2*yBN1 - 2*yBN - 2*h*H*(yBN^2))/h^2 - 2*yAN*(yBN)^2 - d*yBN + Z*yUN;

end
