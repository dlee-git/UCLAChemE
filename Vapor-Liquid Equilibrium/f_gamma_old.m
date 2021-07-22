function[f]=f_gamma(T,x)

% Takes 

%System: Water, Ethanol, Acetone, Acetic Acid
%----------------------------------------------
% Ethanol: 1 part CH3 (1), 1 part CH2 (1), 1 part OH (5)
% Water: 1 part H2O (7)
% Acetone: 1 part CH3CO (9), 1 part CH3 (1)
% Acetic Acid: 1 part CH3 (1), 1 part COOH (20)

% Main group|            | Rk   | Qk
%     1     | CH3:       | .9011| .848
%     1     | CH2:       | .6744| .540
%     5     | OH:        |1.0000| 1.200
%     7     | Water:     | .9200| 1.400
%     9     | CH3CO:     |1.6724| 1.488
%     20    | COOH:      |1.3013| 1.224

% Table H.2 for amk:
% m\k |   1         5         7       9        20
% -------------------------------------------------
% 1   |   0       986.50   1318.0   476.40    663.50
% 5   | 156.40       0     353.50   84.00     199.00        
% 7   | 300.00   -229.10    0      -195.40    -14.090
% 9   | 26.76     164.50   472.50     0       669.40
% 20  | 315.30   -151.00   -66.170 -297.80       0

% UNIFAC Parameters
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

% a_mk = [a11 a15 a17 a19 a120
% a51 a55 a57 a59 a520
% a71 a75 a77 a79 a720
% a91 a95 a97 a99 a920
% a201 a205 a207 a209 a2020]

amk = [0, 986.50, 1318.0, 476.40, 663.50;
156.40, 0, 353.50, 84.00, 199.00;
300.00, -229.10, 0, -195.40, -14.090;
26.76, 164.50, 472.50, 0, 669.40;
315.30, -151.00, -66.170, -297.80, 0];

% calculate ri, qi, and taumk first
ri = [Rk1CH3+Rk1CH2+Rk5, Rk7, Rk9+Rk1CH3, Rk1CH3+Rk20];
qi = [Qk1CH3+Qk1CH2+Qk5, Qk7, Qk9+Qk1CH3, Qk1CH3+Qk20];

%tau_mk = [tau55...tau2020]
tau = exp(-amk/T);

%for J and L
Ji = ri/sum(ri.*x);
Li = qi/sum(qi.*x);

%eki   check if this right
eki = [(1*Qk1CH3+1*Qk1CH2)/qi(1), 0, 1*Qk1CH3/qi(3), 1*Qk1CH3/qi(4);
    1*Qk5/qi(1), 0 , 0 , 0;
    0, Qk7/qi(2), 0 , 0;
    0, 0, Qk9/qi(3), 0;
    0, 0, 0, Qk20/qi(4)];

%beta_ik
beta = zeros(5,5);
for c = 1:5
    for d = 1:5
        beta(c,d) = eki(1,d)*tau(1,c)+eki(2,d)*tau(2,c)+eki(3,d)*tau(3,c)+eki(4,d)*tau(4,c)+eki(5,d)*tau(5,c);
    end
end


%theta_k
theta_k = [(x(1)*qi(1)*eki(1,1)+x(2)*qi(2)*eki(1,2)+x(3)*qi(3)*eki(1,3)+x(4)*qi(4)*eki(1,4))/(x*qi'),...
    (x(1)*qi(1)*eki(2,1)+x(2)*qi(2)*eki(2,2)+x(3)*qi(3)*eki(2,3)+x(4)*qi(4)*eki(2,4))/(x*qi'),...
    (x(1)*qi(1)*eki(3,1)+x(2)*qi(2)*eki(3,2)+x(3)*qi(3)*eki(3,3)+x(4)*qi(4)*eki(3,4))/(x*qi'),...
    (x(1)*qi(1)*eki(4,1)+x(2)*qi(2)*eki(4,2)+x(3)*qi(3)*eki(4,3)+x(4)*qi(4)*eki(4,4))/(x*qi')];

%s_k   // check why its 4x4 but tau is 5x5
s_k = [theta_k(1)*tau(1,1)+theta_k(2)*tau(2,1)+theta_k(3)*tau(3,1)+theta_k(4)*tau(4,1),...
    theta_k(1)*tau(1,2)+theta_k(2)*tau(2,2)+theta_k(3)*tau(3,2)+theta_k(4)*tau(4,2),...
    theta_k(1)*tau(1,3)+theta_k(2)*tau(2,3)+theta_k(3)*tau(3,3)+theta_k(4)*tau(4,3),...
    theta_k(1)*tau(1,4)+theta_k(2)*tau(2,4)+theta_k(3)*tau(3,4)+theta_k(4)*tau(4,4)];

% put everything together!
gC = [exp(1-Ji(1)+log(Ji(1))-5*qi(1)*(1-Ji(1)/Li(1)+log(Ji(1)/Li(1)))),...
    exp(1-Ji(2)+log(Ji(2))-5*qi(2)*(1-Ji(2)/Li(2)+log(Ji(2)/Li(2)))),...
    exp(1-Ji(3)+log(Ji(3))-5*qi(3)*(1-Ji(3)/Li(3)+log(Ji(3)/Li(3)))),...
    exp(1-Ji(4)+log(Ji(4))-5*qi(4)*(1-Ji(4)/Li(4)+log(Ji(4)/Li(4))))];

gR = [exp(qi(1)*(1-...
    (theta_k(1)*beta(1,1)/s_k(1)-eki(1,1)*log(beta(1,1)/s_k(1))+...
    theta_k(2)*beta(2,1)/s_k(2)-eki(2,1)*log(beta(2,1)/s_k(2))+...
    theta_k(3)*beta(3,1)/s_k(3)-eki(3,1)*log(beta(3,1)/s_k(3))+...
    theta_k(4)*beta(4,1)/s_k(4)-eki(4,1)*log(beta(4,1)/s_k(4))))),...
    exp(qi(2)*(1-...
    (theta_k(1)*beta(1,2)/s_k(1)-eki(1,2)*log(beta(1,2)/s_k(1))+...
    theta_k(2)*beta(2,2)/s_k(2)-eki(2,2)*log(beta(2,2)/s_k(2))+...
    theta_k(3)*beta(3,2)/s_k(3)-eki(3,2)*log(beta(3,2)/s_k(3))+...
    theta_k(4)*beta(4,2)/s_k(4)-eki(4,2)*log(beta(4,2)/s_k(4))))),...
    exp(qi(3)*(1-...
    (theta_k(1)*beta(1,3)/s_k(1)-eki(1,3)*log(beta(1,3)/s_k(1))+...
    theta_k(2)*beta(2,3)/s_k(2)-eki(2,3)*log(beta(2,3)/s_k(2))+...
    theta_k(3)*beta(3,3)/s_k(3)-eki(3,3)*log(beta(3,3)/s_k(3))+...
    theta_k(4)*beta(4,3)/s_k(4)-eki(4,3)*log(beta(4,3)/s_k(4))))),...
    exp(qi(4)*(1-...
    (theta_k(1)*beta(1,4)/s_k(1)-eki(1,4)*log(beta(1,4)/s_k(1))+...
    theta_k(2)*beta(2,4)/s_k(2)-eki(2,4)*log(beta(2,4)/s_k(2))+...
    theta_k(3)*beta(3,4)/s_k(3)-eki(3,4)*log(beta(3,4)/s_k(3))+...
    theta_k(4)*beta(4,4)/s_k(4)-eki(4,4)*log(beta(4,4)/s_k(4)))))];
    
    %Send back results
    g = gC.*gR;
    f = g';
    







