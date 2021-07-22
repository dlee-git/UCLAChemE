function [f]=f_Tsat(P)
% f_Psat function calculates the saturated presure as a 4 by 1 matrix given
% a temperature T in kelvin
% Order of rows: Water, Ethanol, Acetone, Acetic Acid
% Order of columns: A, B, C

P = P*101.325;

Antoine = [16.3872, 3885.70, 230.170;
           16.8958, 3795.17, 230.918;
           14.3145, 2756.22, 228.060;
           15.0717, 3580.80, 224.650];

Tsat = (Antoine(:,2)./(Antoine(:,1)-log(P)))-Antoine(:,3) + 273;

f = transpose(Tsat);
end
