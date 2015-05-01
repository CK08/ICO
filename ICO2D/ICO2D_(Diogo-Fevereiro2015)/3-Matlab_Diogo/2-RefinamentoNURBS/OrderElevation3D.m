function [nhx,Uh,nhy,Vh,nhz,Wh,Qw] = OrderElevation3D(n,p,U,m,q,V,f,g,W,Pw,dir,t)

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

bpts = zeros(p+1,q+1,g+1,4);
ebpts = zeros(p+t+1,q+t+1,g+t+1,4);
Nextbpts = zeros(p-1,q-1,g-1,4);
alphas = zeros(p-1,q-1,g-1,1);

nhx = 0;
nhy = 0;
nhz = 0;

if f==1 && g==1 && W==1
    [nc,mc,nds] = size(Pw);
    rc=1;
else
    [nc,mc,rc,nds] = size(Pw);
end

if (dir==1)
    
    bezalfs = zeros(p+t+1,p+1);
    
    n = nc - 1;
    ii = n + p + 1;
    ph = p + t;
    ph2 = floor(ph / 2);
    
    Vh = V;
    Wh = W;
    
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
    
    Qw(1,:,:,:) = Pw(1,:,:,:);
    
    for i=0:ph
        Uh(i+1) = ua;
    end
    
    %Initialize first bezier segment
    for i=0:p
        for row=1:mc
            for row2=1:rc
                bpts(i+1,row,row2,:) = Pw(i+1,row,row2,:);
            end
        end
    end
    
    while(b < ii)
        i = b;
        while b < ii && U(b+1) == U(b+2)
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
                    for row=1:mc
                        for row2=1:rc
                            bpts(k+1,row,row2,:)=alfs(k-s+1)*bpts(k+1,row,row2,:)+(1.0-alfs(k-s+1))*bpts(k,row,row2,:);
                        end
                    end
                end
                
                for row=1:mc
                    for row2=1:rc
                        Nextbpts(saved+1,row,row2,:) = bpts(p+1,row,row2,:);
                    end
                end
            end
        end
        
        %Degree elevation of the Bezier
        for i=lbz:ph
            for row=1:mc
                for row2=1:rc
                    ebpts(i+1,row,row2,:) = 0.0;
                end
            end
            mpi = min(p,i);
            for j=max(0,i-t):mpi
                for row=1:mc
                    for row2=1:rc
                        ebpts(i+1,row,row2,:) = ebpts(i+1,row,row2,:) + bezalfs(i+1,j+1)*bpts(j+1,row,row2,:);
                    end
                end
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
                        for row=1:mc
                            for row2=1:rc
                                Qw(i+1,row,row2,:) = alf*Qw(i+1,row,row2,:) + (1.0-alf)*Qw(i,row,row2,:);
                            end
                        end
                    end
                    
                    if j >=lbz
                        if j-tr <= kind-ph+oldr
                            gam = (ub-Uh(j-tr+1))/den;
                            for row=1:mc
                                for row2=1:rc
                                    ebpts(kj+1,row,row2,:) = gam*ebpts(kj+1,row,row2,:)+(1.0-gam)*ebpts(kj+2,row,row2,:);
                                end
                            end
                        else
                            for row=1:mc
                                for row2=1:rc
                                    ebpts(kj+1,row,row2,:) = bet*ebpts(kj+1,row,row2,:)+(1.0-bet)*ebpts(kj+2,row,row2,:);
                                end
                            end
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
            for row=1:mc
                for row2=1:rc
                    Qw(cind+1,row,row2,:)=ebpts(j+1,row,row2,:);
                end
            end
            cind = cind + 1;
        end
        
        %Set for the next pass
        if b<ii
            for j=0:r-1
                for row=1:mc
                    for row2=1:rc
                        bpts(j+1,row,row2,:)=Nextbpts(j+1,row,row2,:);
                    end
                end
            end
            
            for j=r:p
                for row=1:mc
                    for row2=1:rc
                        bpts(j+1,row,row2,:)=Pw(b-p+j+1,row,row2,:);
                    end
                end
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
    
    nhx = mh-ph-1;


    
%--------------------------------------------------------------------------    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
 
elseif(dir==2)
    
    bezalfs = zeros(q+t+1,q+1);
    
    nn = mc - 1;
    ii = nn + q + 1;
    ph = q + t;
    ph2 = floor(ph / 2);
    
    Uh = U;
    Wh = W;
    
    %Compute Bezier degree elevation coefficients
    bezalfs(1,1) = 1.0;
    bezalfs(ph + 1,q + 1) = 1.0;
    
    for i=1:ph2
        inv = 1.0/Bin(ph,i);
        mpi=min(q,i);
        for j=max(0,i-t):mpi
            bezalfs(i+1,j+1) = inv*Bin(q,j)*Bin(t,i-j);
        end
    end
    
    for i=ph2+1:ph-1
        mpi=min(q,i);
        for j=max(0,i-t):mpi
            bezalfs(i+1,j+1) = bezalfs(ph-i+1,q-j+1);
        end
    end
    
    mh = ph;
    kind = ph + 1;
    r = -1;
    a = q;
    b = q + 1;
    cind = 1;
    
    ua = V(1);
    
    Qw(:,1,:,:) = Pw(:,1,:,:);
    
    for i=0:ph
        Vh(i+1) = ua;
    end
    
    %Initialize first bezier segment
    for i=0:q
        for row=1:nc
            for row2=1:rc
                bpts(row,i+1,row2,:) = Pw(row,i+1,row2,:);
            end
        end
    end
    
    while(b < ii)
        i = b;
        while b < ii && V(b+1) == V(b+2)
            b = b + 1;
        end
        
        mul = b - i + 1;
        mh = mh + mul + t;
        ub = V(b+1);
        oldr = r;
        r = q - mul;
        
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
            for k=q:-1:mul+1
                alfs(k-mul) = numer/(V(a+k+1)-ua);
            end
            
            for j=1:r
                saved = r-j;
                s = mul + j;
                
                for k=q:-1:s
                    for row=1:nc
                        for row2=1:rc
                            bpts(row,k+1,row2,:)=alfs(k-s+1)*bpts(row,k+1,row2,:)+(1.0-alfs(k-s+1))*bpts(row,k,row2,:);
                        end
                    end
                end
                
                for row=1:nc
                    for row2=1:rc
                        Nextbpts(row,saved+1,row2,:) = bpts(row,q+1,row2,:);
                    end
                end
            end
        end
        
        %Degree elevation of the Bezier
        for i=lbz:ph
            for row=1:nc
                for row2=1:rc
                    ebpts(row,i+1,row2,:) = 0.0;
                end
            end
            mpi = min(q,i);
            for j=max(0,i-t):mpi
                for row=1:nc
                    for row2=1:rc
                        ebpts(row,i+1,row2,:) = ebpts(row,i+1,row2,:) + bezalfs(i+1,j+1)*bpts(row,j+1,row2,:);
                    end
                end
            end
        end
        
        
        %Remove inserted knot
        if(oldr > 1)
            first = kind-2;
            last = kind;
            den = ub-ua;
            bet =floor((ub-Vh(kind-1+1))/den);
            
            for tr=1:oldr-1
                i = first;
                j = last;
                kj = j - kind + 1;
                while j-i > tr
                    if(i < cind)
                        alf = (ub-Vh(i+1))/(ua-Vh(i+1));
                        for row=1:nc
                            for row2=1:rc
                                Qw(row,i+1,row2,:) = alf*Qw(row,i+1,row2,:) + (1.0-alf)*Qw(row,i,row2,:);
                            end
                        end
                    end
                    
                    if j >=lbz
                        if j-tr <= kind-ph+oldr
                            gam = (ub-Vh(j-tr+1))/den;
                            for row=1:nc
                                for row2=1:rc
                                    ebpts(row,kj+1,row2,:) = gam*ebpts(row,kj+1,row2,:)+(1.0-gam)*ebpts(row,kj+2,row2,:);
                                end
                            end
                        else
                            for row=1:nc
                                for row2=1:rc
                                    ebpts(row,kj+1,row2,:) = bet*ebpts(row,kj+1,row2,:)+(1.0-bet)*ebpts(row,kj+2,row2,:);
                                end
                            end
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
        if(a ~= q)
            for i=0:ph-oldr-1
                Vh(kind+1) = ua;
                kind = kind + 1;
            end
        end
        
        %Load control points into Qw
        for j=lbz:rbz
            for row=1:nc
                for row2=1:rc
                    Qw(row,cind+1,row2,:)=ebpts(row,j+1,row2,:);
                end
            end
            cind = cind + 1;
        end
        
        %Set for the next pass
        if b<ii
            for j=0:r-1
                for row=1:nc
                    for row2=1:rc
                        bpts(row,j+1,row2,:)=Nextbpts(row,j+1,row2,:);
                    end
                end
            end
            
            for j=r:q
                for row=1:nc
                    for row2=1:rc
                        bpts(row,j+1,row2,:)=Pw(row,b-q+j+1,row2,:);
                    end
                end
            end
            
            a = b;
            b = b + 1;
            ua = ub;
        else
            for i=0:ph
                Vh(kind+i+1)=ub;
            end
        end
        
    end
    
    nhy = mh-ph-1;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


elseif(dir==3)
    
    bezalfs = zeros(g+t+1,g+1);
    
    nn = rc - 1;
    ii = nn + g + 1;
    ph = g + t;
    ph2 = floor(ph / 2);
    
    Uh = U;
    Vh = V;
    
    %Compute Bezier degree elevation coefficients
    bezalfs(1,1) = 1.0;
    bezalfs(ph + 1,g + 1) = 1.0;
    
    for i=1:ph2
        inv = 1.0/Bin(ph,i);
        mpi=min(g,i);
        for j=max(0,i-t):mpi
            bezalfs(i+1,j+1) = inv*Bin(g,j)*Bin(t,i-j);
        end
    end
    
    for i=ph2+1:ph-1
        mpi=min(g,i);
        for j=max(0,i-t):mpi
            bezalfs(i+1,j+1) = bezalfs(ph-i+1,g-j+1);
        end
    end
    
    mh = ph;
    kind = ph + 1;
    r = -1;
    a = g;
    b = g + 1;
    cind = 1;
    
    ua = W(1);
    
    Qw(:,:,1,:) = Pw(:,:,1,:);
    
    for i=0:ph
        Wh(i+1) = ua;
    end
    
    %Initialize first bezier segment
    for i=0:g
        for row=1:nc
            for row2=1:mc
                bpts(row,row2,i+1,:) = Pw(row,row2,i+1,:);
            end
        end
    end
    
    while(b < ii)
        i = b;
        while b < ii && W(b+1) == W(b+2)
            b = b + 1;
        end
        
        mul = b - i + 1;
        mh = mh + mul + t;
        ub = W(b+1);
        oldr = r;
        r = g - mul;
        
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
            for k=g:-1:mul+1
                alfs(k-mul) = numer/(W(a+k+1)-ua);
            end
            
            for j=1:r
                saved = r-j;
                s = mul + j;
                
                for k=g:-1:s
                    for row=1:nc
                        for row2=1:mc
                            bpts(row,row2,k+1,:)=alfs(k-s+1)*bpts(row,row2,k+1,:)+(1.0-alfs(k-s+1))*bpts(row,row2,k,:);
                        end
                    end
                end
                
                for row=1:nc
                    for row2=1:mc
                        Nextbpts(row,row2,saved+1,:) = bpts(row,row2,g+1,:);
                    end
                end
            end
        end
        
        %Degree elevation of the Bezier
        for i=lbz:ph
            for row=1:nc
                for row2=1:mc
                    ebpts(row,row2,i+1,:) = 0.0;
                end
            end
            mpi = min(g,i);
            for j=max(0,i-t):mpi
                for row=1:nc
                    for row2=1:mc
                        ebpts(row,row2,i+1,:) = ebpts(row,row2,i+1,:) + bezalfs(i+1,j+1)*bpts(row,row2,j+1,:);
                    end
                end
            end
        end
        
        
        %Remove inserted knot
        if(oldr > 1)
            first = kind-2;
            last = kind;
            den = ub-ua;
            bet =floor((ub-Wh(kind-1+1))/den);
            
            for tr=1:oldr-1
                i = first;
                j = last;
                kj = j - kind + 1;
                while j-i > tr
                    if(i < cind)
                        alf = (ub-Wh(i+1))/(ua-Wh(i+1));
                        for row=1:nc
                            for row2=1:mc
                                Qw(row,row2,i+1,:) = alf*Qw(row,row2,i+1,:) + (1.0-alf)*Qw(row,row2,i,:);
                            end
                        end
                    end
                    
                    if j >=lbz
                        if j-tr <= kind-ph+oldr
                            gam = (ub-Wh(j-tr+1))/den;
                            for row=1:nc
                                for row2=1:mc
                                    ebpts(row,row2,kj+1,:) = gam*ebpts(row,row2,kj+1,:)+(1.0-gam)*ebpts(row,row2,kj+2,:);
                                end
                            end
                        else
                            for row=1:nc
                                for row2=1:mc
                                    ebpts(row,row2,kj+1,:) = bet*ebpts(row,row2,kj+1,:)+(1.0-bet)*ebpts(row,row2,kj+2,:);
                                end
                            end
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
        if(a ~= g)
            for i=0:ph-oldr-1
                Wh(kind+1) = ua;
                kind = kind + 1;
            end
        end
        
        %Load control points into Qw
        for j=lbz:rbz
            for row=1:nc
                for row2=1:mc
                    Qw(row,row2,cind+1,:)=ebpts(row,row2,j+1,:);
                end
            end
            cind = cind + 1;
        end
        
        %Set for the next pass
        if b<ii
            for j=0:r-1
                for row=1:nc
                    for row2=1:mc
                        bpts(row,row2,j+1,:)=Nextbpts(row,row2,j+1,:);
                    end
                end
            end
            
            for j=r:g
                for row=1:nc
                    for row2=1:mc
                        bpts(row,row2,j+1,:)=Pw(row,row2,b-g+j+1,:);
                    end
                end
            end
            
            a = b;
            b = b + 1;
            ua = ub;
        else
            for i=0:ph
                Wh(kind+i+1)=ub;
            end
        end
        
    end
    
    nhz = mh-ph-1;    

end









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

