clear all
clc

Ns = zeros(12,1);

Ns(1,1) = 0.25;
Ns(3,1) = 0.5;
Ns(5,1) = 0.25;

Ns(7,1) = 0.33333;
Ns(9,1) = 0.33333;
Ns(11,1)= 0.33333;

M = Ns*Ns'