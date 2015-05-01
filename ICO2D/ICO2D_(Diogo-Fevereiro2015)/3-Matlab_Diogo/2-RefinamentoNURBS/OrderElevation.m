function [newStrucure] = OrderElevation (surf, newP, direction)
%--------------------------------------------------------------------------
%   Order elevation of a curve
%   Based on the algorithm A5.9 from "The NURBS book" and functions
%   OrderElevation1D, OrderElevation2D and OrderElevation3D
%   developed by J.F.M. Caseiro at the University of Aveiro.
%
%   Input:
%           - nurbs structure (see nrbmak from the NURBStoolbox) 
%           - order to elevate (integer)
%           - direction (xi:1, eta:2, zeta:3) for the surface and volume
%             cases
%   Output:
%           - new nurbs structure
%
%
%
%   D. Cardoso (diogocardoso@ua.pt)
%   GRIDS
%   University of Aveiro
%--------------------------------------------------------------------------
nargs = nargin;
if nargs < 2
    error ('Not enought arguments while calling the function RefineKnotVec');
elseif nargs>3
    error('Too many input arguments for the function RefineKnotVec');
end

if ~iscell(surf.knots) % is a curve
    % Data from structure
    U = surf.knots;
    Pw= (surf.coefs).';
    n = surf.number;
    p = surf.order-1;
    
    [~,Ubar,Qw] = OrderElevation1D(n,p,U,Pw,newP);
    newStrucure = nrbmak(Qw.', Ubar);
    
elseif iscell(surf.knots) && max(size(surf.knots))==2 % is a surface
    % Data from structure
    U = surf.knots{1};
    V = surf.knots{2};
    Pw= transformPoints(1, surf.coefs);
    n = surf.number(1);
    m = surf.number(2);
    p = surf.order(1)-1;
    q = surf.order(2)-1;
    
    [~,Ubar,~,Vbar,Qw] = OrderElevation2D(n,p,U,m,q,V,Pw,direction,newP);
    Q = transformPoints(3, Qw);
    newStrucure = nrbmak(Q, {Ubar,Vbar});
    
elseif iscell(surf.knots) && max(size(surf.knots))==3 % is a volume
    % Data from structure
    U = surf.knots{1};
    V = surf.knots{2};
    W = surf.knots{3};
    Pw= transformPoints(2, surf.coefs);
    n = surf.number(1);
    m = surf.number(2);
    f = surf.number(3);
    p = surf.order(1)-1;
    q = surf.order(2)-1;
    g = surf.order(3)-1;
    
    [~,Ubar,~,Vbar,~,Wbar,Qw] = OrderElevation3D(n,p,U,m,q,V,f,g,W,Pw,direction,newP);
    Q = transformPoints(4, Qw);
    newStrucure = nrbmak(Q, {Ubar,Vbar,Wbar});
end



end

% Control points transformation matrix
function Pw = transformPoints(type, newPoints)
%
% transform the coefs matrix in list of points
%

% test
% newPoints=surf.coefs;
% newPoints=Qw;[na,nb,~]=size(Qw);

if type==1 %surface
    [a, b, c] = size(newPoints);
    count = 1;
    Pw = zeros(b,c,a);
    for i = 1:c
        for j=1:b
            Pw(j,i,1)=newPoints(1,j,i);
            Pw(j,i,2)=newPoints(2,j,i);
            Pw(j,i,3)=newPoints(3,j,i);
            Pw(j,i,4)=newPoints(4,j,i);
            count = count+1;
        end
    end
    
elseif type==2 % volume
    [a, b, c, d] = size(newPoints);
    count = 1;
    Pw = zeros(b,c,d,a);
    for k=1:d
        for i = 1:c
            for j=1:b
                Pw(j,i,k,1)=newPoints(1,j,i,k);
                Pw(j,i,k,2)=newPoints(2,j,i,k);
                Pw(j,i,k,3)=newPoints(3,j,i,k);
                Pw(j,i,k,4)=newPoints(4,j,i,k);
                count = count+1;
            end
        end
    end
    
elseif type==3 % surface inverse
    
    [na,nb,~]=size(newPoints);
    Pw=zeros(4,na,nb);
    count = 1;
    for i=1:na
        for j=1:nb
            Pw(1,i,j) = newPoints(i,j,1);
            Pw(2,i,j) = newPoints(i,j,2);
            Pw(3,i,j) = newPoints(i,j,3);
            Pw(4,i,j) = newPoints(i,j,4);
            count = count +1;
        end
    end
    
    
elseif type==4 % volume inverse
    
    [na,nb,nc,~]=size(newPoints);
    count = 1;
    Pw=zeros(4,na,nb,nc);
    
    for i=1:na
        for j=1:nb
            for k=1:nc
                Pw(1,i,j,k) = newPoints(i,j,k,1);
                Pw(2,i,j,k) = newPoints(i,j,k,2);
                Pw(3,i,j,k) = newPoints(i,j,k,3);
                Pw(4,i,j,k) = newPoints(i,j,k,4);
                count = count +1;
            end
        end
    end
    
end

end