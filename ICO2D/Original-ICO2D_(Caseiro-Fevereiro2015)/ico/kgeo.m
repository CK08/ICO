clear all
clc

syms g m A k 
syms NNt NDt DNt DDt NTt TNt TTt

syms a11 b11 cm l

%Dimitri

NDt = 1/A*NTt - g/A*NNt;

DNt = 1/A*TNt - g/A*NNt;

DDt = 1/A/A*TTt - g/A/A*TNt - g/A/A*NTt + g*g/A/A*NNt;

matrix = g/m*NNt + (1-g*k/m)*NDt + (1-g*k/m)*DNt + (g*k*k/m-k)*DDt;

%pretty(collect(collect(collect(collect(matrix,NNt),TNt),NTt),TTt));


%Matzen

A = cm;
k = b11;
m = l*l;

NDt = 1/A*NTt - g/A*NNt;

DNt = 1/A*TNt - g/A*NNt;

DDt = 1/A/A*TTt - g/A/A*TNt - g/A/A*NTt + g*g/A/A*NNt;

matrix = g/m*NNt + (1-g*k/m)*NDt + (1-g*k/m)*DNt + (g*k*k/m-k)*DDt;

pretty(collect(collect(collect(collect(matrix,NNt),TNt),NTt),TTt));