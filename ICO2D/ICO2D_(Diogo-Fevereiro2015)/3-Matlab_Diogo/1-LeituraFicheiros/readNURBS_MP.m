function [ncp, orders, knots, coefs] = readNURBS_MP(fid)
%----------------------------------------------------------------
%
% Read function for the NURBSInput.txt and NURBSOutput.txt file 
% from ICO2D
%
% Generates four variables which are:
%       - ncp (Number of control points per direction)
%       - orders (Order in each direction)
%       - knots (knot vectors in each direction)
%       - coefs (the control points)
% each structure/variable contains two fields for each patch
%
% This feature will be updated in future versions so that the 
% number of patches is created dynamically.
%
% Copyright (C) Diogo Cardoso, 2015
%
%----------------------------------------------------------------

% Read the header

npatch=fgetl(fid);
ncpf1 = str2num(fgetl(fid));
ncpf2 = str2num(fgetl(fid));
ordersf1 = str2num(fgetl(fid));
ordersf2 = str2num(fgetl(fid));

if (max(size(ordersf1))==3) && (max(size(ordersf2))==3)
    % Patch 1
    Uf1 = str2num(fgetl(fid));
    Vf1 = str2num(fgetl(fid));
    Wf1 = str2num(fgetl(fid));
    knotsf1 = {Uf1,Vf1,Wf1};
    coefsf1 = zeros(4,ncpf1(1),ncpf1(2),ncpf1(3));
    
    for i = 1:ncpf1(1)
        for j = 1:ncpf1(2)
            for k = 1:ncpf1(3)
                coor = str2num(fgetl(fid));
                coefsf1(1,i,j,k) = coor(1)*coor(4);
                coefsf1(2,i,j,k) = coor(2)*coor(4);
                coefsf1(3,i,j,k) = coor(3)*coor(4);
                coefsf1(4,i,j,k) = coor(4);
            end
        end
    end
    % Patch 2
    Uf2 = str2num(fgetl(fid));
    Vf2 = str2num(fgetl(fid));
    Wf2 = str2num(fgetl(fid));
    knotsf2 = {Uf2,Vf2,Wf2};
    coefsf2 = zeros(4,ncpf2(1),ncpf2(2),ncpf2(3));
    
    for i = 1:ncpf2(1)
        for j = 1:ncpf2(2)
            for k = 1:ncpf2(3)
                coor = str2num(fgetl(fid));
                coefsf2(1,i,j,k) = coor(1)*coor(4);
                coefsf2(2,i,j,k) = coor(2)*coor(4);
                coefsf2(3,i,j,k) = coor(3)*coor(4);
                coefsf2(4,i,j,k) = coor(4);
            end
        end
    end
    
elseif (max(size(ordersf1))==2) && (max(size(ordersf2))==2)
    % Patch 1
    Uf1 = str2num(fgetl(fid));
    Vf1 = str2num(fgetl(fid));
    knotsf1 = {Uf1,Vf1};
    coefsf1 = zeros(3,ncpf1(1),ncpf1(2));
    
    for i = 1:ncpf1(1)
        for j = 1:ncpf1(2)
            
            coor = str2num(fgetl(fid));
            coefsf1(1,i,j) = coor(1)*coor(3);
            coefsf1(2,i,j) = coor(2)*coor(3);
            coefsf1(3,i,j) = coor(3);
            
        end
    end
    % Patch 2
    Uf2 = str2num(fgetl(fid));
    Vf2 = str2num(fgetl(fid));
    knotsf2 = {Uf2,Vf2};
    coefsf2 = zeros(3,ncpf2(1),ncpf2(2));
    
    for i = 1:ncpf2(1)
        for j = 1:ncpf2(2)
            
            coor = str2num(fgetl(fid));
            coefsf2(1,i,j) = coor(1)*coor(3);
            coefsf2(2,i,j) = coor(2)*coor(3);
            coefsf2(3,i,j) = coor(3);
            
        end
    end
end

% Aloc to output variables
ncp{1} = ncpf1; ncp{2} = ncpf2;
orders{1} = ordersf1; orders{2} = ordersf2;
knots{1} = knotsf1; knots{2} = knotsf2;
coefs{1} = coefsf1; coefs{2} = coefsf2;

end