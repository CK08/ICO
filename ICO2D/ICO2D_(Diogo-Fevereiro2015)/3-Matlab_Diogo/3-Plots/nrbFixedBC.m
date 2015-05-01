function [x,y,z] = nrbFixedBC(nurbs,data,scale)
%
% Representation of a nurbs surface or volume with the corresponding
% Dirichlet boundary conditions (homogeneous).
%
% Input:
%       - nurbs structure (see nrbmak);
%       - data: [n x 2] matrix with the direction of the boundary condition
%       to apply in the second column and the corresponding index point on
%       the first;
%       - scale: [x,y,z] vector with the scale of the pyramid to represent
%       the boundary condition.
%
% Output: 
%       - [x,y,z] for posterior fill() operations;
%
% Copyright (C) Diogo Cardoso 2014
%

coor = nrbind2coor(nurbs,data(:,1));


if max(size(nurbs.order))==3
    np = [100 100 100];
elseif max(size(nurbs.order))==2
    np = [100 100];
else
    return;
end

nrbplotx(nurbs,np);hold on;
[x, y, z] = BC(coor,scale, data(:,2));
fill3(x,y,z,'b');
view([1 1 1]);axis square;


end



function [x,y,z]=BC(points,scale,direction)
%
% Points must be a matrix with n-lines of points and three columns (x,y,z)
% Scale of arrow:  vecto with three values [scale X , scale Y, scale Z] -
% default [1, 1, 1]
% Directions must be a matrix for each point, defined as X = 1 , Y = 2 and Z = 3
%           or X = -1, Y = -2 and Z = -3 for negative directions

n=size(points,1);
ni=1;
nf=5;
scaleX=scale(1);
scaleY=scale(2);
scaleZ=scale(3);
% For each point
for p=1:n
    ponta=points(p,:);
    
    % Ox positivo
    if direction(p,1)==1
        zz=[ponta(3)-scaleZ/2 ponta(3)+scaleZ/2 ponta(3)+scaleZ/2 ponta(3)-scaleZ/2;
            ponta(3)-scaleZ/2 ponta(3)+scaleZ/2 ponta(3) ponta(3);
            ponta(3)+scaleZ/2 ponta(3)+scaleZ/2 ponta(3) ponta(3);
            ponta(3)+scaleZ/2 ponta(3)-scaleZ/2 ponta(3) ponta(3);
            ponta(3)-scaleZ/2 ponta(3)-scaleZ/2 ponta(3) ponta(3)];
        zz=zz.';
        yy=[ponta(2)-scaleY/2 ponta(2)-scaleY/2 ponta(2)+scaleY/2 ponta(2)+scaleY/2;
            ponta(2)-scaleY/2 ponta(2)-scaleY/2 ponta(2) ponta(2);
            ponta(2)-scaleY/2 ponta(2)+scaleY/2 ponta(2) ponta(2);
            ponta(2)+scaleY/2 ponta(2)+scaleY/2 ponta(2) ponta(2);
            ponta(2)+scaleY/2 ponta(2)-scaleY/2 ponta(2) ponta(2)];
        yy=yy.';
        xx=[ponta(1)-scaleX ponta(1)-scaleX ponta(1)-scaleX ponta(1)-scaleX;
            ponta(1)-scaleX ponta(1)-scaleX ponta(1) ponta(1);
            ponta(1)-scaleX ponta(1)-scaleX ponta(1) ponta(1);
            ponta(1)-scaleX ponta(1)-scaleX ponta(1) ponta(1);
            ponta(1)-scaleX ponta(1)-scaleX ponta(1) ponta(1)];
        xx=xx.';
    % Ox negativo
    elseif direction(p,1)==-1
        zz=[ponta(3)-scaleZ/2 ponta(3)+scaleZ/2 ponta(3)+scaleZ/2 ponta(3)-scaleZ/2;
            ponta(3)-scaleZ/2 ponta(3)+scaleZ/2 ponta(3) ponta(3);
            ponta(3)+scaleZ/2 ponta(3)+scaleZ/2 ponta(3) ponta(3);
            ponta(3)+scaleZ/2 ponta(3)-scaleZ/2 ponta(3) ponta(3);
            ponta(3)-scaleZ/2 ponta(3)-scaleZ/2 ponta(3) ponta(3)];
        zz=zz.';
        yy=[ponta(2)-scaleY/2 ponta(2)-scaleY/2 ponta(2)+scaleY/2 ponta(2)+scaleY/2;
            ponta(2)-scaleY/2 ponta(2)-scaleY/2 ponta(2) ponta(2);
            ponta(2)-scaleY/2 ponta(2)+scaleY/2 ponta(2) ponta(2);
            ponta(2)+scaleY/2 ponta(2)+scaleY/2 ponta(2) ponta(2);
            ponta(2)+scaleY/2 ponta(2)-scaleY/2 ponta(2) ponta(2)];
        yy=yy.';
        scaleX=-scaleX;
        xx=[ponta(1)-scaleX ponta(1)-scaleX ponta(1)-scaleX ponta(1)-scaleX;
            ponta(1)-scaleX ponta(1)-scaleX ponta(1) ponta(1);
            ponta(1)-scaleX ponta(1)-scaleX ponta(1) ponta(1);
            ponta(1)-scaleX ponta(1)-scaleX ponta(1) ponta(1);
            ponta(1)-scaleX ponta(1)-scaleX ponta(1) ponta(1)];
        xx=xx.';    
    % Oy positivo
    elseif direction(p,1)==2
        zz=[ponta(3)-scaleZ/2 ponta(3)+scaleZ/2 ponta(3)+scaleZ/2 ponta(3)-scaleZ/2;
            ponta(3)-scaleZ/2 ponta(3)+scaleZ/2 ponta(3) ponta(3);
            ponta(3)+scaleZ/2 ponta(3)+scaleZ/2 ponta(3) ponta(3);
            ponta(3)+scaleZ/2 ponta(3)-scaleZ/2 ponta(3) ponta(3);
            ponta(3)-scaleZ/2 ponta(3)-scaleZ/2 ponta(3) ponta(3)];
        zz=zz.';
        xx=[ponta(1)-scaleX/2 ponta(1)-scaleX/2 ponta(1)+scaleX/2 ponta(1)+scaleX/2;
            ponta(1)-scaleX/2 ponta(1)-scaleX/2 ponta(1) ponta(1);
            ponta(1)-scaleX/2 ponta(1)+scaleX/2 ponta(1) ponta(1);
            ponta(1)+scaleX/2 ponta(1)+scaleX/2 ponta(1) ponta(1);
            ponta(1)+scaleX/2 ponta(1)-scaleX/2 ponta(1) ponta(1)];
        xx=xx.';
        yy=[ponta(2)-scaleY ponta(2)-scaleY ponta(2)-scaleY ponta(2)-scaleY;
            ponta(2)-scaleY ponta(2)-scaleY ponta(2) ponta(2);
            ponta(2)-scaleY ponta(2)-scaleY ponta(2) ponta(2);
            ponta(2)-scaleY ponta(2)-scaleY ponta(2) ponta(2);
            ponta(2)-scaleY ponta(2)-scaleY ponta(2) ponta(2)];
        yy=yy.';
    % Oy negativo
    elseif direction(p,1)==-2
        zz=[ponta(3)-scaleZ/2 ponta(3)+scaleZ/2 ponta(3)+scaleZ/2 ponta(3)-scaleZ/2;
            ponta(3)-scaleZ/2 ponta(3)+scaleZ/2 ponta(3) ponta(3);
            ponta(3)+scaleZ/2 ponta(3)+scaleZ/2 ponta(3) ponta(3);
            ponta(3)+scaleZ/2 ponta(3)-scaleZ/2 ponta(3) ponta(3);
            ponta(3)-scaleZ/2 ponta(3)-scaleZ/2 ponta(3) ponta(3)];
        zz=zz.';
        xx=[ponta(1)-scaleX/2 ponta(1)-scaleX/2 ponta(1)+scaleX/2 ponta(1)+scaleX/2;
            ponta(1)-scaleX/2 ponta(1)-scaleX/2 ponta(1) ponta(1);
            ponta(1)-scaleX/2 ponta(1)+scaleX/2 ponta(1) ponta(1);
            ponta(1)+scaleX/2 ponta(1)+scaleX/2 ponta(1) ponta(1);
            ponta(1)+scaleX/2 ponta(1)-scaleX/2 ponta(1) ponta(1)];
        xx=xx.';
        scaleY=-scaleY;
        yy=[ponta(2)-scaleY ponta(2)-scaleY ponta(2)-scaleY ponta(2)-scaleY;
            ponta(2)-scaleY ponta(2)-scaleY ponta(2) ponta(2);
            ponta(2)-scaleY ponta(2)-scaleY ponta(2) ponta(2);
            ponta(2)-scaleY ponta(2)-scaleY ponta(2) ponta(2);
            ponta(2)-scaleY ponta(2)-scaleY ponta(2) ponta(2)];
        yy=yy.';
    % Oz positivo
    elseif direction(p,1)==3
        xx=[ponta(1)-scaleX/2 ponta(1)+scaleX/2 ponta(1)+scaleX/2 ponta(1)-scaleX/2;
            ponta(1)-scaleX/2 ponta(1)+scaleX/2 ponta(1) ponta(1);
            ponta(1)+scaleX/2 ponta(1)+scaleX/2 ponta(1) ponta(1);
            ponta(1)+scaleX/2 ponta(1)-scaleX/2 ponta(1) ponta(1);
            ponta(1)-scaleX/2 ponta(1)-scaleX/2 ponta(1) ponta(1)];
        xx=xx.';
        yy=[ponta(2)-scaleY/2 ponta(2)-scaleY/2 ponta(2)+scaleY/2 ponta(2)+scaleY/2;
            ponta(2)-scaleY/2 ponta(2)-scaleY/2 ponta(2) ponta(2);
            ponta(2)-scaleY/2 ponta(2)+scaleY/2 ponta(2) ponta(2);
            ponta(2)+scaleY/2 ponta(2)+scaleY/2 ponta(2) ponta(2);
            ponta(2)+scaleY/2 ponta(2)-scaleY/2 ponta(2) ponta(2)];
        yy=yy.';
        zz=[ponta(3)-scaleZ ponta(3)-scaleZ ponta(3)-scaleZ ponta(3)-scaleZ;
            ponta(3)-scaleZ ponta(3)-scaleZ ponta(3) ponta(3);
            ponta(3)-scaleZ ponta(3)-scaleZ ponta(3) ponta(3);
            ponta(3)-scaleZ ponta(3)-scaleZ ponta(3) ponta(3);
            ponta(3)-scaleZ ponta(3)-scaleZ ponta(3) ponta(3)];
        zz=zz.';
    % Oz negativo
    elseif direction(p,1)==-3
        xx=[ponta(1)-scaleX/2 ponta(1)+scaleX/2 ponta(1)+scaleX/2 ponta(1)-scaleX/2;
            ponta(1)-scaleX/2 ponta(1)+scaleX/2 ponta(1) ponta(1);
            ponta(1)+scaleX/2 ponta(1)+scaleX/2 ponta(1) ponta(1);
            ponta(1)+scaleX/2 ponta(1)-scaleX/2 ponta(1) ponta(1);
            ponta(1)-scaleX/2 ponta(1)-scaleX/2 ponta(1) ponta(1)];
        xx=xx.';
        yy=[ponta(2)-scaleY/2 ponta(2)-scaleY/2 ponta(2)+scaleY/2 ponta(2)+scaleY/2;
            ponta(2)-scaleY/2 ponta(2)-scaleY/2 ponta(2) ponta(2);
            ponta(2)-scaleY/2 ponta(2)+scaleY/2 ponta(2) ponta(2);
            ponta(2)+scaleY/2 ponta(2)+scaleY/2 ponta(2) ponta(2);
            ponta(2)+scaleY/2 ponta(2)-scaleY/2 ponta(2) ponta(2)];
        yy=yy.';
        scaleZ=-scaleZ;
        zz=[ponta(3)-scaleZ ponta(3)-scaleZ ponta(3)-scaleZ ponta(3)-scaleZ;
            ponta(3)-scaleZ ponta(3)-scaleZ ponta(3) ponta(3);
            ponta(3)-scaleZ ponta(3)-scaleZ ponta(3) ponta(3);
            ponta(3)-scaleZ ponta(3)-scaleZ ponta(3) ponta(3);
            ponta(3)-scaleZ ponta(3)-scaleZ ponta(3) ponta(3)];
        zz=zz.';
    end
    x(:,ni:nf)=xx;
    y(:,ni:nf)=yy;
    z(:,ni:nf)=zz;
    ni=nf+1;
    nf=nf+5;
end
end