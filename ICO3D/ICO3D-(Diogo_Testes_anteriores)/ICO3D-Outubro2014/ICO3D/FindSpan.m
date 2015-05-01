function [s] = FindSpan(n,p,u,U)

%--------------------------------------------------------------------------
%   Calculate knot span to which the data point u belongs
%   Based on the algorithm A2.1 from "The NURBS book"
%    
%   Input:
%       n = nknots - p - 1
%       p - B-Spline degree
%       u - data point
%       U - knot vector
%   Output:
%       s - Knot span of data point u
%--------------------------------------------------------------------------

%Special Case - The point is the last knot in the knot span
if(u == U(n+1)) 
    s=n;
    return
end

low = p;
high = n+1;
mid = low + 1;

%high = n+1;
%mid = (low + high)/2;

while(u < U(mid) || u >= U(mid+1))
   if(u < U(mid)) 
       mid = mid - 1;
   else
       low = mid;
       mid = mid + 1;
   end
end
    
s = mid;
return