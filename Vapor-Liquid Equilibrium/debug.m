% Takes 
T = 298.15;
x = [0.25, 0.25, 0.25, 0.25];


Rk1CH3 = .9011;
Rk1CH2 = .6744;
Rk5 = 1.0000;
Rk7 = .9200;
Rk9 = 1.6724;
Rk20 = 1.3013;
Qk1CH3 = .848;
Qk1CH2 = .540;
Qk5 = 1.200;
Qk7 = 1.400;
Qk9 = 1.488;
Qk20 = 1.224;

amk = [0, 0, 986.50, 1318.0, 476.40, 663.50;
0, 0, 986.50, 1318.0, 476.40, 663.50;
156.40, 156.4, 0, 353.50, 84.00, 199.00;
300.00, 300.00, -229.10, 0, -195.40, -14.090;
26.76, 26.76, 164.50, 472.50, 0, 669.40;
315.30, 315.30, -151.00, -66.170, -297.80, 0];

% calculate ri, qi, and taumk first
ri = [Rk1CH3+Rk1CH2+Rk5, Rk7, Rk9+Rk1CH3, Rk1CH3+Rk20];
qi = [Qk1CH3+Qk1CH2+Qk5, Qk7, Qk9+Qk1CH3, Qk1CH3+Qk20];

%tau_mk = [tau55...tau2020]
tau = exp(-amk/T);

%for J and L
Ji = ri./sum(ri.*x);
Li = qi./sum(qi.*x);

eki = [(Qk1CH3)/qi(1), 0, 1*Qk1CH3/qi(3), 1*Qk1CH3/qi(4);
    Qk1CH2/qi(1), 0, 0, 0;
    Qk5/qi(1), 0 , 0 , 0;
    0, Qk7/qi(2), 0 , 0;
    0, 0, Qk9/qi(3), 0;
    0, 0, 0, Qk20/qi(4)];

%beta_ik
beta = zeros(6,4);
for subgroup = 1:6
    for species = 1:4
        beta(subgroup,species) = eki(1,species)*tau(1,subgroup)+eki(2,species)*tau(2,subgroup)+eki(3,species)*tau(3,subgroup)+eki(4,species)*tau(4,subgroup)+eki(5,species)*tau(5,subgroup)+eki(6,species)*tau(6,subgroup);
    end
end

theta_k = [(x(1)*qi(1)*eki(1,1)+x(2)*qi(2)*eki(1,2)+x(3)*qi(3)*eki(1,3)+x(4)*qi(4)*eki(1,4))/sum(sum((x.*qi))),...
    (x(1)*qi(1)*eki(2,1)+x(2)*qi(2)*eki(2,2)+x(3)*qi(3)*eki(2,3)+x(4)*qi(4)*eki(2,4))/sum((x.*qi)),...
    (x(1)*qi(1)*eki(3,1)+x(2)*qi(2)*eki(3,2)+x(3)*qi(3)*eki(3,3)+x(4)*qi(4)*eki(3,4))/sum((x.*qi)),...
    (x(1)*qi(1)*eki(4,1)+x(2)*qi(2)*eki(4,2)+x(3)*qi(3)*eki(4,3)+x(4)*qi(4)*eki(4,4))/sum((x.*qi)),...
    (x(1)*qi(1)*eki(5,1)+x(2)*qi(2)*eki(5,2)+x(3)*qi(3)*eki(5,3)+x(4)*qi(4)*eki(5,4))/sum((x.*qi)),...
    (x(1)*qi(1)*eki(6,1)+x(2)*qi(2)*eki(6,2)+x(3)*qi(3)*eki(6,3)+x(4)*qi(4)*eki(6,4))/sum((x.*qi))];

s_k = [theta_k(1)*tau(1,1)+theta_k(2)*tau(2,1)+theta_k(3)*tau(3,1)+theta_k(4)*tau(4,1)+theta_k(5)*tau(5,1)+theta_k(6)*tau(6,1),...
       theta_k(1)*tau(1,2)+theta_k(2)*tau(2,2)+theta_k(3)*tau(3,2)+theta_k(4)*tau(4,2)+theta_k(5)*tau(5,2)+theta_k(6)*tau(6,2),...
       theta_k(1)*tau(1,3)+theta_k(2)*tau(2,3)+theta_k(3)*tau(3,3)+theta_k(4)*tau(4,3)+theta_k(5)*tau(5,3)+theta_k(6)*tau(6,3),...
       theta_k(1)*tau(1,4)+theta_k(2)*tau(2,4)+theta_k(3)*tau(3,4)+theta_k(4)*tau(4,4)+theta_k(5)*tau(5,4)+theta_k(6)*tau(6,4),...
       theta_k(1)*tau(1,5)+theta_k(2)*tau(2,5)+theta_k(3)*tau(3,5)+theta_k(4)*tau(4,5)+theta_k(5)*tau(5,5)+theta_k(6)*tau(6,5),...
       theta_k(1)*tau(1,6)+theta_k(2)*tau(2,6)+theta_k(3)*tau(3,6)+theta_k(4)*tau(4,6)+theta_k(5)*tau(5,6)+theta_k(6)*tau(6,6)];

gC = [exp(1-Ji(1)+log(Ji(1))-5*qi(1)*(1-Ji(1)/Li(1)+log(Ji(1)/Li(1)))),...
    exp(1-Ji(2)+log(Ji(2))-5*qi(2)*(1-Ji(2)/Li(2)+log(Ji(2)/Li(2)))),...
    exp(1-Ji(3)+log(Ji(3))-5*qi(3)*(1-Ji(3)/Li(3)+log(Ji(3)/Li(3)))),...
    exp(1-Ji(4)+log(Ji(4))-5*qi(4)*(1-Ji(4)/Li(4)+log(Ji(4)/Li(4))))];

gR = [exp(qi(1)*(1-...
    (theta_k(1)*beta(1,1)/s_k(1)-eki(1,1)*log(beta(1,1)/s_k(1))+...
    theta_k(2)*beta(2,1)/s_k(2)-eki(2,1)*log(beta(2,1)/s_k(2))+...
    theta_k(3)*beta(3,1)/s_k(3)-eki(3,1)*log(beta(3,1)/s_k(3))+...
    theta_k(4)*beta(4,1)/s_k(4)-eki(4,1)*log(beta(4,1)/s_k(4))+...
    theta_k(5)*beta(5,1)/s_k(5)-eki(5,1)*log(beta(5,1)/s_k(5))+...
    theta_k(6)*beta(6,1)/s_k(6)-eki(6,1)*log(beta(6,1)/s_k(6))))),...
    exp(qi(2)*(1-...
    (theta_k(1)*beta(1,2)/s_k(1)-eki(1,2)*log(beta(1,2)/s_k(1))+...
    theta_k(2)*beta(2,2)/s_k(2)-eki(2,2)*log(beta(2,2)/s_k(2))+...
    theta_k(3)*beta(3,2)/s_k(3)-eki(3,2)*log(beta(3,2)/s_k(3))+...
    theta_k(4)*beta(4,2)/s_k(4)-eki(4,2)*log(beta(4,2)/s_k(4))+...
    theta_k(5)*beta(5,2)/s_k(5)-eki(5,2)*log(beta(5,2)/s_k(5))+...
    theta_k(6)*beta(6,2)/s_k(6)-eki(6,2)*log(beta(6,2)/s_k(6))))),...
    exp(qi(3)*(1-...
    (theta_k(1)*beta(1,3)/s_k(1)-eki(1,3)*log(beta(1,3)/s_k(1))+...
    theta_k(2)*beta(2,3)/s_k(2)-eki(2,3)*log(beta(2,3)/s_k(2))+...
    theta_k(3)*beta(3,3)/s_k(3)-eki(3,3)*log(beta(3,3)/s_k(3))+...
    theta_k(4)*beta(4,3)/s_k(4)-eki(4,3)*log(beta(4,3)/s_k(4))+...
    theta_k(5)*beta(5,3)/s_k(5)-eki(5,3)*log(beta(5,3)/s_k(5))+...
    theta_k(6)*beta(6,3)/s_k(6)-eki(6,3)*log(beta(6,3)/s_k(6))))),...
    exp(qi(4)*(1-...
    (theta_k(1)*beta(1,4)/s_k(1)-eki(1,4)*log(beta(1,4)/s_k(1))+...
    theta_k(2)*beta(2,4)/s_k(2)-eki(2,4)*log(beta(2,4)/s_k(2))+...
    theta_k(3)*beta(3,4)/s_k(3)-eki(3,4)*log(beta(3,4)/s_k(3))+...
    theta_k(4)*beta(4,4)/s_k(4)-eki(4,4)*log(beta(4,4)/s_k(4))+...
    theta_k(5)*beta(5,4)/s_k(5)-eki(5,4)*log(beta(5,4)/s_k(5))+...
    theta_k(6)*beta(6,4)/s_k(6)-eki(6,4)*log(beta(6,4)/s_k(6)))))];
 

    %Send back results
    g = gC.*gR;
