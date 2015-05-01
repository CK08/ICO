function [nh,Uh,Qw] = OrderElevation1D(n,p,U,Pw,t)

%--------------------------------------------------------------------------
%   Order elevation of a curve
%   Based on the algorithm A5.9 from "The NURBS book"
%    
%   Input:
%       n = number of basis functions
%       p - B-Spline degree
%       U - knot vector
%       Pw - Control points pre-multiplied by the weight
%       t - degree to elevate
%   Output:
%       nh - new number of basis functions
%       Uh - New knot vector
%       Qw - New control points multiplied by new weight
%
%   J.F.M. Caseiro
%   GRIDS
%   University of Aveiro
%--------------------------------------------------------------------------

bezalfs = zeros(p+t+1,p+1);
bpts = zeros(p+1,3);
ebpts = zeros(p+t+1,3);
Nextbpts = zeros(p-1,3);
alphas = zeros(p-1,1);

[nc,mc] = size(Pw);

n = nc - 1;

m = n + p + 1;
ph = p + t;
ph2 = floor(ph / 2);


%Compute Bezier degree elevation coefficients
bezalfs(1,1) = 1.0;
bezalfs(ph + 1,p + 1) = 1.0;

for i=1:ph2
    inv = 1.0/Bin(ph,i);
    mpi=min(p,i);
    for j=max(0,i-t):mpi
        bezalfs(i+1,j+1) = inv*Bin(p,j)*Bin(t,i-j);
    end
end

for i=ph2+1:ph-1
    mpi=min(p,i);
    for j=max(0,i-t):mpi
        bezalfs(i+1,j+1) = bezalfs(ph-i+1,p-j+1);
    end
end

mh = ph;
kind = ph + 1;
r = -1;
a = p;
b = p + 1;
cind = 1;

ua = U(1);

Qw(1,:) = Pw(1,:);

for i=0:ph
    Uh(i+1) = ua;
end

%Initialize first bezier segment
for i=0:p
    bpts(i+1,:) = Pw(i+1,:);
end

while(b < m)
    i = b;
    while b < m && U(b+1) == U(b+2)
        b = b + 1;
    end
    
    mul = b - i + 1;
    mh = mh + mul + t;
    ub = U(b+1);
    oldr = r;
    r = p - mul;
    
    %Insert knot u(b) r times
    if(oldr > 0)
        lbz = floor((oldr+2)/2);
    else
        lbz = 1;
    end
    
    if(r > 0)
        rbz = ph-floor((r+1)/2);
    else
        rbz = ph;
    end
    
    if(r > 0)
        numer = ub-ua;
        for k=p:-1:mul+1
            alfs(k-mul) = numer/(U(a+k+1)-ua);
        end
        
        for j=1:r
            saved = r-j;
            s = mul + j;
            
            for k=p:-1:s
                bpts(k+1,:)=alfs(k-s+1)*bpts(k+1,:)+(1.0-alfs(k-s+1))*bpts(k,:);
            end
            Nextbpts(saved+1,:) = bpts(p+1,:);
        end
    end
    
    %Degree elevation of the Bezier
    for i=lbz:ph
        ebpts(i+1,:) = 0.0;
        mpi = min(p,i);
        for j=max(0,i-t):mpi
            ebpts(i+1,:) = ebpts(i+1,:) + bezalfs(i+1,j+1)*bpts(j+1,:);
        end
    end
    
    
    %Remove inserted knot
    if(oldr > 1)
        first = kind-2;
        last = kind;
        den = ub-ua;
        bet =floor((ub-Uh(kind-1+1))/den);
        
        for tr=1:oldr-1
            i = first;
            j = last;
            kj = j - kind + 1;
            while j-i > tr
                if(i < cind)
                    alf = (ub-Uh(i+1))/(ua-Uh(i+1));
                    Qw(i+1,:) = alf*Qw(i+1,:) + (1.0-alf)*Qw(i,:);
                end
                
                if j >=lbz
                    if j-tr <= kind-ph+oldr
                        gam = (ub-Uh(j-tr+1))/den;
                        ebpts(kj+1,:) = gam*ebpts(kj+1,:)+(1.0-gam)*ebpts(kj+2,:);
                    else
                        ebpts(kj+1,:) = bet*ebpts(kj+1,:)+(1.0-bet)*ebpts(kj+2,:);
                    end
                end
                
                i = i + 1;
                j = j - 1;
                kj = kj - 1;
                
            end
            
            first = first - 1;
            last = last + 1;
            
        end
    end
    
    %Load the knot ua
    if(a ~= p)
        for i=0:ph-oldr-1
            Uh(kind+1) = ua;
            kind = kind + 1;
        end
    end
    
    %Load control points into Qw
    for j=lbz:rbz
        Qw(cind+1,:)=ebpts(j+1,:);
        cind = cind + 1;
    end
    
    %Set for the next pass
    if b<m
        for j=0:r-1
            bpts(j+1,:)=Nextbpts(j+1,:);
        end
        
        for j=r:p
            bpts(j+1,:)=Pw(b-p+j+1,:);
        end
        
        a = b;
        b = b + 1;
        ua = ub;
    else
        for i=0:ph
            Uh(kind+i+1)=ub;
        end
    end
    
    
end


nh = mh-ph-1;










end


function b = Bin(n,k)
%  Computes the binomial coefficient.
%
%      ( n )      n!
%      (   ) = --------
%      ( k )   k!(n-k)!
%
%  b = bincoeff(n,k)
%
%  Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215.

                                                          % double bincoeff(int n, int k)
                                                          % {
b = floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));      %   return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
end 

function f = factln(n)
% computes ln(n!)
if n <= 1, f = 0; return, end
f = gammaln(n+1); %log(factorial(n));
end

