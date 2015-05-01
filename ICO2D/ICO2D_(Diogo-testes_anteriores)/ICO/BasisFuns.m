function [BF] = BasisFuns(s,u,p,U)

%--------------------------------------------------------------------------
%   Calculate basis functions of a B-Spline
%   Based on the algorithm A2.2 from "The NURBS book"
%    
%   Input:
%       s - knot span
%       u - data point
%       p - B-Spline degree
%       U - knot vector
%   Output:
%       BF(p+1,1) - Non-zero basis functions
%--------------------------------------------------------------------------


N = zeros(p+1);
BF = zeros(p+1,1);
left = zeros(p+1);
right = zeros(p+1);

N(1) = 1.0;

i = s;

for j=1:p
    left(j+1) = u - U(i+1-j);
    right(j+1) = U(i+j) - u;
    saved = 0.0;
    
    for r=0:j-1
        temp = N(r+1)/(right(r+2)+left(j-r+1));
        N(r+1) = saved+right(r+2)*temp;
        saved = left(j-r+1)*temp;
    end
    N(j+1) = saved;
    continue
end

for k1=1:p+1
    BF(k1,1) = N(k1);
end

return
end
