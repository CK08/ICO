function [Ubar,Qw] = RefineKnotVecCurve (n,p,U,Pw,X,r)

%--------------------------------------------------------------------------
%   Knot refinement of a curve
%   Based on the algorithm A5.4 from "The NURBS book"
%    
%   Input:
%       n = nknots - p - 1
%       p - B-Spline degree
%       U - knot vector
%       Pw - Control points pre-multiplied by the weight
%       X - knot vector to be inserted
%       r - length of X
%   Output:
%       Ubar - New knot vector
%       Qw - New control points multiplied by new weight 
%
%   J.F.M. Caseiro
%   GRIDS
%   University of Aveiro
%--------------------------------------------------------------------------

n = n - 1;
r = r - 1;

m = n + p + 1;

a = FindSpan(n,p,X(1),U);
b = FindSpan(n,p,X(r+1),U);

%%%%%%%
a = a - 1;
%b = b + 1;
%%%%%%

for j=0:a-p
    Qw(j+1,:) = Pw(j+1,:);
end

for j=b-1:n
    Qw(j+r+2,:) = Pw(j+1,:);
end

for j=0:a
    Ubar(j+1) = U(j+1);
end

for j=b+p:m
    Ubar(j+r+2) = U(j+1);
end

i = b + p - 1;
k = b + p + r;

for j=r:-1:0
    
    while(X(j+1) <= U(i+1) && i>a)
        Qw(k-p-1+1,:) = Pw(i-p-1+1,:);
        Ubar(k+1) = U(i+1);
        k = k-1;
        i = i-1;
    end
    
    Qw(k-p-1+1,:) = Qw(k-p+1,:);
    
    for l=1:p
        
        ind = k-p+l;
        alfa = Ubar(k+l+1)-X(j+1);
        if(abs(alfa)==0.0)
            Qw(ind-1+1,:) = Qw(ind+1,:);
        else
            alfa = alfa/(Ubar(k+l+1)-U(i-p+l+1));
            Qw(ind-1+1,:) = alfa*Qw(ind-1+1,:)+(1.0-alfa)*Qw(ind+1,:);
        end
    end
    
    Ubar(k+1) = X(j+1);
    k = k-1;
    
end


end