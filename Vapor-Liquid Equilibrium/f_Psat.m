function [f]=f_Psat(T)
% f_Psat function calculates the saturated presure as a 4 by 1 matrix given
% a temperature T in kelvin
% Order of rows: Water, Ethanol, Acetone, Acetic Acid
% Order of columns: A, B, C

T = T - 273;

% Antoine Constants
Antoine = [16.3872, 3885.70, 230.170;
           16.8958, 3795.17, 230.918;
           14.3145, 2756.22, 228.060;
           15.0717, 3580.80, 224.650];

% Calculate Psat
Psat = exp(Antoine(:,1)-Antoine(:,2)./(T+Antoine(:,3)))/101.325;

% Correctly orient output
f = transpose(Psat);
end