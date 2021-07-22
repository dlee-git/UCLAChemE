clear all; close all; clc;

% Given parameters
T = 273;

x = [40, 24, 13, 23];

% Initial assumption(s): phi=1
phi= [1;1];

% Generating initial guess for pressure
psat = f_Psat(T);