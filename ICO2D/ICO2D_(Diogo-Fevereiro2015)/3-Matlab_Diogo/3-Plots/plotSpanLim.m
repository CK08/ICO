function plotSpanLim(p,q,U,V,np,CP,col)
% The plotSpanLim plots the span limits of a NURBS or Bspline patch surface
% Input:
%       - degrees (p,q); 
%       - knot vectors (U,V); 
%       - np points along each edge.
%       - Control Points (coefs): see nrbmak
%       - color (RGB)
% Output:
%       - plot of the limits
%
%
% Copyright(C): Diogo Cardoso, 2014
%
%

% Security check
nargs = nargin;
if nargs == 2 
    if ~isstruct(p)
        error ('Need a NURBS to plot and the number of subdivisions!');
    else
        nurbs = p;
        np = q;
        p = nurbs.order(1)-1;
        q = nurbs.order(2)-1;
        U = nurbs.knots{1};
        V = nurbs.knots{2};
        CP = nurbs.coefs;
        color = 'k';
    end
elseif nargs == 3 
    if ~isstruct(p)
        error ('Need a NURBS to plot and the number of subdivisions!');
    else
        nurbs = p;
        np = q;
        color = U;
        if ~isnumeric(color) && length(color)>1, color = 'k';end
        p = nurbs.order(1)-1;
        q = nurbs.order(2)-1;
        U = nurbs.knots{1};
        V = nurbs.knots{2};
        CP = nurbs.coefs;
        
    end    
elseif nargs>2 && nargs<6
  error ('Parameters are inconsistent')
elseif nargs == 6
    color = 'k';
elseif nargs == 7
    color = col;
end
%-----------

lengU = max(length(U)); 
lengV = max(length(V)); 

nb = lengU-p-1;
nq = lengV-q-1;
    
u = linspace(U(1),U(end),np);
v = linspace(V(1),V(end),np);

[a, b, c] = size(CP);
Pw = zeros(b,c,a);
for i = 1:c
    for j=1:b
        Pw(j,i,1)=CP(1,j,i);
        Pw(j,i,2)=CP(2,j,i);
        Pw(j,i,3)=CP(3,j,i);
        Pw(j,i,4)=CP(4,j,i);
    end
end

SBFu = zeros(p+1,np);

for i=p+1:lengU-p
        
    
        k1=1;
        SRF = zeros(1,np,3);
        
        %Find span index (xi)(starting from 1)
        uspan = FindSpan(nb,p,U(i),U);
        
        %Find basis functions values (xi)
        BFu = BasisFuns(uspan,U(i),p,U);
        SBFu(:,k1) = BFu(:,:);
        
        SBFv = zeros(q+1,np);
        
        for k2=1:np
            %Find span index (eta) (starting from 1)
            vspan = FindSpan(nq,q,v(k2),V);
            
            %Find basis functions values (eta)
            BFv = BasisFuns(vspan,v(k2),q,V);
            SBFv(:,k2) = BFv(:,:);
            
            den = 0.0;
            for l=0:q
                temp = zeros(1,1,4);
                for m=0:p
                    temp  = temp  + BFu(m+1,1)*Pw(uspan-p+m, vspan-q+l, :);
                    continue
                end
                
                den = den + temp(1,1,4)*BFv(l+1,1);
                SRF(k1,k2,1) = SRF(k1,k2,1) + temp(1,1,1)*BFv(l+1,1);
                SRF(k1,k2,2) = SRF(k1,k2,2) + temp(1,1,2)*BFv(l+1,1);
                SRF(k1,k2,3) = SRF(k1,k2,3) + temp(1,1,3)*BFv(l+1,1);
            end
            
            SRF(k1,k2,1) = SRF(k1,k2,1)/den;
            SRF(k1,k2,2) = SRF(k1,k2,2)/den;
            SRF(k1,k2,3) = SRF(k1,k2,3)/den;
            
            
        end
        
        %subplot(spdim,1,umesh);
        hold on;
        plot3(SRF(1,:,1),SRF(1,:,2),SRF(1,:,3),'Color',color,'LineStyle','-','LineWidth',2);        

    
end

for i=q+1:lengV-q
    
        
        k2=1;
        SRF = zeros(np,1,3);
        for k1=1:np
            
        
            %Find span index (xi)(starting from 1)
            uspan = FindSpan(nb,p,u(k1),U);

            %Find basis functions values (xi)
            BFu = BasisFuns(uspan,u(k1),p,U);
            SBFu(:,k1) = BFu(:,:);
        
            
            %for k2=1:np
            %Find span index (eta) (starting from 1)
            vspan = FindSpan(nq,q,V(i),V);
            
            %Find basis functions values (eta)
            BFv = BasisFuns(vspan,V(i),q,V);
            SBFv(:,k2) = BFv(:,:);
            
            den = 0.0;
            for l=0:q
                temp = zeros(1,1,4);
                for m=0:p
                    temp  = temp  + BFu(m+1,1)*Pw(uspan-p+m, vspan-q+l, :);
                    continue
                end
                
                den = den + temp(1,1,4)*BFv(l+1,1);
                SRF(k1,k2,1) = SRF(k1,k2,1) + temp(1,1,1)*BFv(l+1,1);
                SRF(k1,k2,2) = SRF(k1,k2,2) + temp(1,1,2)*BFv(l+1,1);
                SRF(k1,k2,3) = SRF(k1,k2,3) + temp(1,1,3)*BFv(l+1,1);
            end
            
            SRF(k1,k2,1) = SRF(k1,k2,1)/den;
            SRF(k1,k2,2) = SRF(k1,k2,2)/den;
            SRF(k1,k2,3) = SRF(k1,k2,3)/den;
            
        end
        hold on;
        plot3(SRF(:,1,1),SRF(:,1,2),SRF(:,1,3),'Color',color,'LineStyle','-','LineWidth',2);

   
end


end