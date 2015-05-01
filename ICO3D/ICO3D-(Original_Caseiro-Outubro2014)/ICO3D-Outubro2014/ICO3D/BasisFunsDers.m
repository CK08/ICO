function [ders] = BasisFunsDers(s,u,p,n,U)

%--------------------------------------------------------------------------
%   Calculate derivatives of the basis functions of a B-Spline
%   Based on the algorithm A2.3 from "The NURBS book"
%    
%   Input:
%       s - knot span
%       u - data point
%       p - B-Spline degree
%       n - nth derivative
%       U - knot vector
%   Output:
%       ders(n+1,p+1) - Derivative of the basis functions
%                       Line i contains the ith derivative
%--------------------------------------------------------------------------

ders = zeros(n+1,p+1);
ndu = zeros(p+1,p+1);
left = zeros(p+1);
right = zeros(p+1);
a = zeros(2,p+1);

ndu(1,1) = 1.0;

for j=1:p
    left(j+1) = u - U(s+1-j);
    right(j+1) = U(s+j) - u;
    saved = 0.0;
    
    for r=0:j-1
        ndu(j+1,r+1) = right(r+2)+left(j-r+1);
        temp = ndu(r+1,j)/ndu(j+1,r+1);
        
        ndu(r+1,j+1) = saved + right(r+2)*temp;
        saved = left(j-r+1)*temp;
%         temp = N(r+1)/(right(r+2)+left(j-r+1));
%         N(r+1) = saved+right(r+2)*temp;
%         saved = left(j-r+1)*temp;
    end
    ndu(j+1,j+1) = saved;
    continue
end

%Load basis functions
for j=0:p
    ders(1,j+1) = ndu(j+1,p+1);
end

for r=0:p
   s1 = 0;
   s2 = 1;
   a(1,1) = 1.0;
   
   for k=1:n
       d =0.0;
       rk = r-k ;
       pk = p-k;
       
       if(r>=k)
           a(s2+1,1) = a(s1+1,1)/ndu(pk+2,rk+1);
           d = a(s2+1,1)*ndu(rk+1,pk+1);
       end
       
       if(rk>= -1)
           j1 = 1;
       else
           j1 = -rk;
       end
       
       if((r-1) <= pk)
           j2 = k-1;
       else
           j2 = p-r;
       end
      
       for j=j1:j2
           a(s2+1,j+1) = (a(s1+1,j+1)-a(s1+1,j))/ndu(pk+2,rk+j+1);
           d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1);
       end
        
       if(r<=pk)
           a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1);
           d = d + a(s2+1,k+1)*ndu(r+1,pk+1);
       end
       
       ders(k+1,r+1) = d;
       j = s1;
       s1 = s2;
       s2 = j + 1;
        
        
   end
   
end

r = p;
for k = 1:n
    for j = 0:p
        ders(k+1,j+1) = ders(k+1,j+1)*r;
    end
    r = r*(p-k);
end

return