P = 101.325;
T = 298.15;
y = [0.25,0.25,0.25,0.25];

% Order: Water, Ethanol, Acetone, Acetic Acid
% k (empirical interaction parameter specific to an i-j
% molecular pair and for chemically similar species, kij = 0

omega = [0.345, 0.645, 0.307, 0.467];
Tc = [647.1, 513.9, 508.2, 592.0];
Pc = [22055, 6148, 4701, 5786];
Zc = [0.229, 0.240, 0.233, 0.211];
Vc = [55.9, 167, 209, 179.7];
P_sat= [50, 145, 70, 120];

%mixing rules
R = 8314; %[cc-kpa/mol-K)
t = T;
n = length(y);            
omegaij = zeros(n,n);
Tcij = zeros(n,n);
Pcij = zeros(n,n);
Zcij = zeros(n,n);
Vcij = zeros(n,n);
Trij = zeros(n,n);
B0ij = zeros(n,n);
B1ij = zeros(n,n);
Bhij = zeros(n,n);
Bij = zeros(n,n);

% First, use Prausnitz's Mixing Rules
for i = 1:n
    for j = 1:n
        omegaij(i,j) = (omega(i)+omega(j))/2;
        Tcij(i,j) = (Tc(i)*Tc(j))^0.5;
        Zcij(i,j) = (Zc(i)+Zc(j))/2;
        Vcij(i,j) = ((Vc(i)^(1/3)+Vc(j)^(1/3))/2)^3;
        Pcij(i,j) = Zcij(i,j)*R*Tcij(i,j)/Vcij(i,j);
        Trij(i,j) = t/Tcij(i,j);
        B0ij(i,j) = 0.083-0.422/Trij(i,j)^1.6;
        B1ij(i,j) = 0.139-0.172/Trij(i,j)^4.2;
        Bhij(i,j) = B0ij(i,j)+omegaij(i,j)*B1ij(i,j);
        Bij(i,j) = Bhij(i,j)*R*Tcij(i,j)/Pcij(i,j);
    end
end

%Resetting and working using equation 11.64 from the book
B = Bij;    

% n is the number of components, and del is del_ij
n = length(B(1,:));
del = zeros(size(B));

% set values for del_ij
for i = 1:n
    for j = 1:n
        del(i,j) = 2*B(i,j)-B(i,i)-B(j,j);
    end
end

phi = zeros(1,n);

% now evaluate eqn 11.64
for k=1:n
    inner = B(k,k)*(P-P_sat(k));
    for i = 1:n
        for j = 1:n
            inner = inner + 0.5*y(i)*y(j)*(2*del(i,k)-del(i,j))*P;
        end
    end
    phi(k) = exp(inner/(R*T));
end

f = phi';