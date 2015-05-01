clear all
close all
clc

x = rand(9,1)*4-2;
y = rand(9,1)*4-2;
z = x.*exp(-x.^2-y.^2);

F = TriScatteredInterp(x,y,z);
ti = -2:.25:2;
[qx,qy] = meshgrid(ti,ti);
qz = F(qx,qy);

figure(1);scatter3(x,y,z);
hold on;
mesh(qx,qy,qz);

figure(2);
contour(qx,qy,qz);


clear all
close all
clc

x = linspace(-2*pi,2*pi);
y = linspace(0,4*pi);
[X,Y] = meshgrid(x,y);
Z = sin(X)+cos(Y);