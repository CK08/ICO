function [Ubar,Vbar,Wbar,Qw] = RefineKnotVecVolume (n,p,U,m,q,V,f,g,W,Pw,X,r,dir)

%--------------------------------------------------------------------------
%   Knot refinement of a volume
%   Based on the algorithm A5.4 from "The NURBS book"
%    
%   Input:
%       n = nknots_x - p - 1
%       p - B-Spline degree in xx direction
%       U - xi knot vector
%       m = nknots_y - q - 1
%       q - B-Spline degree in yy direction
%       V - eta knot vector
%       f = nknots_z - g - 1
%       g - B-Spline degree in zz direction
%       W - zeta knot vector
%       Pw - Control points pre-multiplied by the weight
%       X - knot vector to be inserted
%       r - length of X
%       dir - direction of the refinement (1-xi, 2-eta, 3-zeta)
%   Output:
%       Ubar - New xi knot vector
%       Vbar - New eta knot vector
%       Wbar - New zeta knot vector
%       Qw - New control points multiplied by new weight
%
%   J.F.M. Caseiro
%   GRIDS
%   University of Aveiro
%--------------------------------------------------------------------------

if(dir == 1)
    
    n = n - 1;
    r = r - 1;
    
    a = FindSpan(n,p,X(1),U);
    b = FindSpan(n,p,X(r+1),U);
    
    a = a - 1;

    Vbar = V;
    Wbar = W;
    
    for j=0:a-p
        Qw(j+1,:,:,:) = Pw(j+1,:,:,:);
    end
    
    for j=b-1:n
        Qw(j+r+2,:,:,:) = Pw(j+1,:,:,:);
    end
    
    for j=0:a
        Ubar(j+1) = U(j+1);
    end
    
    for j=b+p:n + p + 1
        Ubar(j+r+2) = U(j+1);
    end
    
    i = b + p - 1;
    k = b + p + r;
    
    for j=r:-1:0
        
        while(X(j+1) <= U(i+1) && i>a)
            Qw(k-p-1+1,:,:,:) = Pw(i-p-1+1,:,:,:);
            Ubar(k+1) = U(i+1);
            k = k-1;
            i = i-1;
        end
        
        Qw(k-p-1+1,:,:,:) = Qw(k-p+1,:,:,:);
        
        for l=1:p
            
            ind = k-p+l;
            alfa = Ubar(k+l+1)-X(j+1);
            if(abs(alfa)==0.0)
                Qw(ind-1+1,:,:,:) = Qw(ind+1,:,:,:);
            else
                alfa = alfa/(Ubar(k+l+1)-U(i-p+l+1));
                Qw(ind-1+1,:,:,:) = alfa*Qw(ind-1+1,:,:,:)+(1.0-alfa)*Qw(ind+1,:,:,:);
            end
        end
        
        Ubar(k+1) = X(j+1);
        k = k-1;
        
    end

%----------------------------------------------------------------------------------    
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------
    
elseif dir==2
    
    m = m - 1;
    r = r - 1;
    
    a = FindSpan(m,q,X(1),V);
    b = FindSpan(m,q,X(r+1),V);
    
    a = a - 1;

    Ubar = U;
    Wbar = W;
    
    for j=0:a-q
        Qw(:,j+1,:,:) = Pw(:,j+1,:,:);
    end
    
    for j=b-1:m
        Qw(:,j+r+2,:,:) = Pw(:,j+1,:,:);
    end
    
    for j=0:a
        Vbar(j+1) = V(j+1);
    end
    
    for j=b+q:m + q + 1
        Vbar(j+r+2) = V(j+1);
    end
    
    i = b + q - 1;
    k = b + q + r;
    
    for j=r:-1:0
        
        while(X(j+1) <= V(i+1) && i>a)
            Qw(:,k-q-1+1,:,:) = Pw(:,i-q-1+1,:,:);
            Vbar(k+1) = V(i+1);
            k = k-1;
            i = i-1;
        end
        
        Qw(:,k-q-1+1,:,:) = Qw(:,k-q+1,:,:);
        
        for l=1:q
            
            ind = k-q+l;
            alfa = Vbar(k+l+1)-X(j+1);
            if(abs(alfa)==0.0)
                Qw(:,ind-1+1,:,:) = Qw(:,ind+1,:,:);
            else
                alfa = alfa/(Vbar(k+l+1)-V(i-q+l+1));
                Qw(:,ind-1+1,:,:) = alfa*Qw(:,ind-1+1,:,:)+(1.0-alfa)*Qw(:,ind+1,:,:);
            end
        end
        
        Vbar(k+1) = X(j+1);
        k = k-1;
        
    end
    
    
%----------------------------------------------------------------------------------    
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------
%----------------------------------------------------------------------------------
    
elseif dir==3
    
    f = f - 1;
    r = r - 1;
    
    a = FindSpan(f,g,X(1),W);
    b = FindSpan(f,g,X(r+1),W);
    
    a = a - 1;

    Ubar = U;
    Vbar = V;
    
    for j=0:a-g
        Qw(:,:,j+1,:) = Pw(:,:,j+1,:);
    end
    
    for j=b-1:f
        Qw(:,:,j+r+2,:) = Pw(:,:,j+1,:);
    end
    
    for j=0:a
        Wbar(j+1) = W(j+1);
    end
    
    for j=b+g:f + g + 1
        Wbar(j+r+2) = W(j+1);
    end
    
    i = b + g - 1;
    k = b + g + r;
    
    for j=r:-1:0
        
        while(X(j+1) <= W(i+1) && i>a)
            Qw(:,:,k-g-1+1,:) = Pw(:,:,i-g-1+1,:);
            Wbar(k+1) = W(i+1);
            k = k-1;
            i = i-1;
        end
        
        Qw(:,:,k-g-1+1,:) = Qw(:,:,k-g+1,:);
        
        for l=1:g
            
            ind = k-g+l;
            alfa = Wbar(k+l+1)-X(j+1);
            if(abs(alfa)==0.0)
                Qw(:,:,ind-1+1,:) = Qw(:,:,ind+1,:);
            else
                alfa = alfa/(Wbar(k+l+1)-W(i-g+l+1));
                Qw(:,:,ind-1+1,:) = alfa*Qw(:,:,ind-1+1,:)+(1.0-alfa)*Qw(:,:,ind+1,:);
            end
        end
        
        Wbar(k+1) = X(j+1);
        k = k-1;
        
    end
    
    
    
    
end

end


