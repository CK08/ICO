function [newStrucure] = RefineKnotVec(structure, knot, direction)
%--------------------------------------------------------------------------
%   Knot refinement of a curve, surface and volume
%   Based on the algorithm A5.4 from "The NURBS book" and functions
%   RefineKnotVecCurve, RefineKnotVecSurface and RefineKnotVecVolume
%   developed by J.F.M. Caseiro at the University of Aveiro.
%
%   Input:
%           - nurbs structure (see nrbmak from the NURBStoolbox) 
%           - coordinate of the knot or vector of knots coordinates to be 
%             inserted
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

if ~iscell(structure.knots) % is a curve
    % Data from structure
    U = structure.knots;
    Pw= (structure.coefs).';
    r = length(knot);
    n = structure.number;
    p = structure.order-1;
    
    [Ubar,Qw] = RefineKnotVecCurve (n,p,U,Pw,knot,r);
    newStrucure = nrbmak(Qw.', Ubar);
    
elseif iscell(structure.knots) && max(size(structure.knots))==2 % is a surface
    % Data from structure
    U = structure.knots{1};
    V = structure.knots{2};
    Pw= transformPoints(1, structure.coefs);
    r = length(knot);
    n = structure.number(1);
    m = structure.number(2);
    p = structure.order(1)-1;
    q = structure.order(2)-1;
    
    [Ubar,Vbar,Qw] = RefineKnotVecSurface (n,p,U,m,q,V,Pw,knot,r,direction);
    Q = transformPoints(3, {Qw,length(Ubar)-p-1,length(Vbar)-q-1});
    newStrucure = nrbmak(Q, {Ubar,Vbar});
    
elseif iscell(structure.knots) && max(size(structure.knots))==3 % is a volume
    % Data from structure
    U = structure.knots{1};
    V = structure.knots{2};
    W = structure.knots{3};
    Pw= transformPoints(2, structure.coefs);
    r = length(knot);
    n = structure.number(1);
    m = structure.number(2);
    f = structure.number(3);
    p = structure.order(1)-1;
    q = structure.order(2)-1;
    g = structure.order(3)-1;
    
    [Ubar,Vbar,Wbar,Qw] = RefineKnotVecVolume (n,p,U,m,q,V,f,g,W,Pw,knot,r,direction);
    Q = transformPoints(4, {Qw,length(Ubar)-p-1,length(Vbar)-q-1, length(Wbar)-g-1});
    newStrucure = nrbmak(Q, {Ubar,Vbar,Wbar});
end




end

% Control points transformation matrix
function newPoints = transformPoints(type, points)
%
% transform the coefs matrix in list of points
%

if type==1 %surface
    [a, b, c] = size(points);
    count = 1;
    newPoints = zeros(b,c,a);
    for i = 1:c
        for j=1:b
            newPoints(j,i,1)=points(1,j,i);
            newPoints(j,i,2)=points(2,j,i);
            newPoints(j,i,3)=points(3,j,i);
            newPoints(j,i,4)=points(4,j,i);
            count = count+1;
        end
    end
    
elseif type==2 % volume
    [a, b, c, d] = size(points);
    count = 1;
    newPoints = zeros(b,c,d,a);
    for k=1:d
        for i = 1:c
            for j=1:b
                newPoints(j,i,k,1)=points(1,j,i,k);
                newPoints(j,i,k,2)=points(2,j,i,k);
                newPoints(j,i,k,3)=points(3,j,i,k);
                newPoints(j,i,k,4)=points(4,j,i,k);
                count = count+1;
            end
        end
    end
    
elseif type==3 % surface inverse
    Qw = points{1};
    na  = points{2};
    nb  = points{3};
    newPoints=zeros(4,na,nb);
    count = 1;
    for i=1:na
        for j=1:nb
            newPoints(1,i,j) = Qw(i,j,1);
            newPoints(2,i,j) = Qw(i,j,2);
            newPoints(3,i,j) = Qw(i,j,3);
            newPoints(4,i,j) = Qw(i,j,4);
            count = count +1;
        end
    end
    
    
elseif type==4 % volume inverse
    Qw = points{1};
    na  = points{2};
    nb  = points{3};
    nc  = points{4};
    newPoints=zeros(4,na,nb,nc);
    count = 1;
    for i=1:na
        for j=1:nb
            for k=1:nc
                newPoints(1,i,j,k) = Qw(i,j,k,1);
                newPoints(2,i,j,k) = Qw(i,j,k,2);
                newPoints(3,i,j,k) = Qw(i,j,k,3);
                newPoints(4,i,j,k) = Qw(i,j,k,4);
                count = count +1;
            end
        end
    end
    
end

end