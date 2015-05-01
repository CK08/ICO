clear all
clc

x = [-0.57 -0.57 0.57 0.57];
y = [-0.57 0.57 -0.57 0.57];

data = [100.0 25.0 175.0 200.0];

% zgrid = gridfit(x,y,data,4,4)
% [X,Y] = meshgrid(x,y)

% scatter(x,y,50,data,'fill','s')
% colorbar

[xq,yq] = meshgrid(x,y);
vq = griddata(x,y,data,xq,yq)
contour(xq,yq,vq)
colorbar