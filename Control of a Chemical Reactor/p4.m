clc; clear all; close all;
 
Gp = tf([4000 5075.2],[1 538.5193 -37748.894 -37299.32]);
freq = logspace(0, 1, 2500);
bode(Gp, freq)[magnitude, phase, freq] = bode(Gp, freq);

