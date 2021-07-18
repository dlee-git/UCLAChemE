clear all; clc;

%Define Step Size 
nodes = 100; h = 1/(nodes-1);
% Define Constants
d = 0.0; Z = 0.0; D = 0.1; B = 1.5; Y = 0.05; E = 0.0; H = 0.00; theta = 0.0;
% Define Errror
error = 1.0;

% Intialize matrixes for Jacobian and F
J_matrix = zeros(3*nodes-1, 3*nodes-1);
F_vector = zeros(3*nodes-1, 1);
% Define Y which will host all 3 concentrations in an array.
Y_current = zeros(3*nodes-1,1);
Y_next = zeros(3*nodes-1,1);
% Define X - distance from r = 0 to r = 1 discreted into 100 points
X = linspace(0,1,nodes);


%Iterate multivariable Newton's Method until it reaches appropriate error
% Y_M+1 = Y_M - J(Y_M)^-1 * F(Y_M)
while (error > 0.1) 
    %Define Jacobian and F for each iteration 
    %yA
    J_matrix(1,1) = -2*D*(1+h)/h^2 - Y_current(nodes+1,1)^2 - Y;
    J_matrix(1,2) = 2*D/h^2;
    J_matrix(nodes+1,1) = Y_current(nodes+1,1)^2;
    F_vector(1,1) = A_first_function(Y_current(1,1),Y_current(2,1),Y_current(nodes+1,1));
    for i = 2:nodes-1 
        J_matrix(i,i-1)= D/h^2;
        J_matrix(i,i) = -2*D/h^2 - Y_current(i+nodes,1)^2 - Y;
        J_matrix(i,i+1) = D/h^2;
        J_matrix(i+nodes, i) = -2*Y_current(i+nodes,1)^2;
        F_vector(i,1) = A_middle_function(Y_current(i+1,1),Y_current(i,1),Y_current(i-1,1),Y_current(i+nodes,1));
    end
    J_matrix(nodes, nodes-1) = 2*D/h^2;
    J_matrix(nodes, nodes) = -2*D*(1+E*h)/h^2 - Y_current(2*nodes,1)^2 - Y;
    J_matrix(2*nodes, nodes) = -2*Y_current(2*nodes,1)^2;
    F_vector(nodes, 1) = A_last_function(Y_current(nodes-1,1),Y_current(nodes,1), Y_current(2*nodes,1));
    %yB
    J_matrix(1,nodes+1) = -2*Y_current(1,1)*Y_current(nodes+1,1);
    J_matrix(nodes+1, nodes+1) = -2*D*(1+h)/h^2 - 4*Y_current(1,1)*Y_current(nodes+1,1) - d;
    J_matrix(nodes+1,nodes+2) = 2*D/h^2;
    F_vector(nodes+1,1) = B_first_function(Y_current(1,1),Y_current(nodes+1,1),Y_current(nodes+2,1),0);
    for i = (nodes+2):(2*nodes-1)
        J_matrix(i-nodes,i) = -2*Y_current(i-nodes,1)*Y_current(i,1);
        J_matrix(i,i-1) = D/h^2;
        J_matrix(i,i) = -2*D/h^2 - 4*Y_current(i-nodes,1)*Y_current(i,1) - d;
        J_matrix(i,i+1) = D/h^2;
        J_matrix(i+nodes-1,i) = d;
        F_vector(i,1) = B_middle_function(Y_current(i-nodes,1),Y_current(i+1,1),Y_current(i,1),Y_current(i-1,1),Y_current(i+nodes,1));
    end
    J_matrix(nodes,2*nodes) = 2*Y_current(nodes,1)*Y_current(2*nodes,1);
    J_matrix(2*nodes,2*nodes-1) = 2*D/h^2;
    J_matrix(2*nodes,2*nodes) = D*(-2-4*h*H*Y_current(2*nodes,1))/h^2 - 4*Y_current(nodes,1)*Y_current(2*nodes) - d;
    J_matrix(3*nodes-1,2*nodes) = d;
    F_vector(2*nodes,1) = B_last_function(Y_current(nodes,1),Y_current(2*nodes-1,1),Y_current(2*nodes,1), Y_current(3*nodes-1,1));
    %yU
    J_matrix(nodes+2,2*nodes+1) = Z;
    J_matrix(2*nodes+1,2*nodes+1) = -2*D/h^2 - Z;
    J_matrix(2*nodes+1,2*nodes+2) = D/h^2;
    F_vector(2*nodes+1,1) = U_middle_function(Y_current(nodes+1,1),Y_current(2*nodes+2,1),Y_current(2*nodes+1,1),0);
    for i = (2*nodes+2):(3*nodes-2)
        J_matrix(i-nodes+1,i) = Z;
        J_matrix(i,i-1) = D/h^2;
        J_matrix(i,i) = -2*D/h^2 - Z;
        J_matrix(i,i+1) = D/h^2;
        F_vector(i,1) = U_middle_function(Y_current(i-nodes,1),Y_current(i+1,1),Y_current(i,1),Y_current(i-1,1));
    end
    J_matrix(2*nodes,3*nodes-1) = Z;
    J_matrix(3*nodes-1,3*nodes-2) = 2*D/h^2;
    J_matrix(3*nodes-1,3*nodes-1)= -2*D*(1+h*theta)/h^2 - Z;
    F_vector(3*nodes-1,1) = U_last_function(Y_current(2*nodes,1),Y_current(3*nodes-2,1),Y_current(3*nodes-1,1));
    %Calculate new ym1
    Y_next = Y_current - inv(J_matrix)*F_vector;
    %Determine if we should keep iterating 
    error = 0;
    for j = 1:3*nodes-1
       if F_vector(j,1)^2 > 0.0001 
           error = 1;
       end
    end
    %Update values
    Y_current = Y_next;
end
%Create matrices to hold concentration values for each species
yA(:,1) = Y_current(1:nodes,1);
yB(:,1) = Y_current(nodes+1:2*nodes,1);
yU(2:nodes,1) = Y_current(2*nodes+1:3*nodes-1,1);
yU(1,1) = 0;
%Plot Graphs
figure(1);
plot(X, yA, 'g');
hold on
plot(X, yB, 'r');
plot(X, yU, 'b');
xlabel('x');
ylabel('y');
title('Steady-State Concentration Profiles Q.A');
legend('yA', 'yB', 'yU');

% Transient State with Explict Euler Method
delta_t = 0.0001;

% Define Y which will host all 3 concentrations in an array.
Yt_current = zeros(3*nodes-1,1);
Yt_next = zeros(3*nodes-1,1);

err = 1;
count = 0;

while err > 0.01
    %Use Explict Euler to determine new y_transient values
    for i = 1:3*nodes-1 
       if i == 1
           Yt_next(1,1) = Yt_current(1,1) + delta_t*A_first_function(Yt_current(1,1),Yt_current(2,1),Yt_current(nodes+1,1));
       elseif i < nodes
           Yt_next(i,1) = Yt_current(i,1) + delta_t*A_middle_function(Yt_current(i+1,1),Yt_current(i,1),Yt_current(i-1,1),Yt_current(i+nodes,1));
       elseif i == nodes
           Yt_next(i,1) = Yt_current(i,1) + delta_t*A_last_function(Yt_current(i-1,1),Yt_current(i,1),Yt_current(2*i,1));
       elseif i == nodes+1 
           Yt_next(i,1) = Yt_current(i,1) + delta_t*B_first_function(Yt_current(1,1),Yt_current(i,1),Yt_current(i+1,1),0);
       elseif i < 2*nodes
           Yt_next(i,1) = Yt_current(i,1) + delta_t*B_middle_function(Yt_current(i-nodes,1),Yt_current(i+1,1),Yt_current(i,1),Yt_current(i-1,1),Yt_current(i+nodes,1));
       elseif i == 2*nodes
           Yt_next(i,1) = Yt_current(i,1) + delta_t*B_last_function(Yt_current(nodes,1),Yt_current(i-1,1),Yt_current(i,1), Yt_current(3*nodes-1,1));
       elseif i == 2*nodes+1 
           Yt_next(i,1) = Yt_current(i,1) + delta_t*U_middle_function(Yt_current(i-nodes,1),Yt_current(i+1,1),Yt_current(i,1),0);
       elseif i < 3*nodes-1
           Yt_next(i,1) = Yt_current(i,1) + delta_t*U_middle_function(Yt_current(i-nodes,1),Yt_current(i+1,1),Yt_current(i,1),Yt_current(i-1,1));
       elseif i == 3*nodes-1
           Yt_next(i,1) = Yt_current(i,1) + delta_t*U_last_function(Yt_current(2*nodes,1),Yt_current(i-1,1),Yt_current(i,1));
       end
    end
    %Determine error for this iteration 
    err = norm(abs(Yt_next-Y_current))/norm(Y_current);
    %Plot curves
    count = count + 1;
    % Plot every 20000 iterations
    if mod(count,20000) == 0
        yAt(:,1) = Yt_next(1:nodes,1);
        yBt(:,1) = Yt_next(nodes+1:2*nodes,1);
        yUt(1,1) = 0;
        yUt(2:nodes,1) = Yt_next(2*nodes+1:3*nodes-1,1);
        
        figure(2);
        plot(X,yAt);
        hold on;
        figure(3);
        plot(X,yBt);
        hold on;
        figure(4);
        yUt(1,1) = 0;
        plot(X,yUt);
        hold on;
    end
    Yt_current = Yt_next;
end
t = delta_t*count; 

figure(2);
xlabel('x');
ylabel('yA');
title('Evolution of Spatial Profiles of yA');

figure(3);
xlabel('x');
ylabel('yB');
title('Evolution of Spatial Profiles of yB');

figure(4);
xlabel('x');
ylabel('yU');
title('Evolution of Spatial Profiles of yU');

fprintf('Time elapsed to convergence: %f \n', t);
