%Lagrange Multiplier Example

clear all
clc

%Minimise the function
%

%Subject to the restriction g>=0
%g = u-1;

u = 1;
lm = -1;

for i=1:10
    
    g = u-1;
    f = 1/4*u^4;
    
    if(g<= 0 || lm <=0)
    
        Kc = [3*u^2, 1; 1, 0];
        Fc = [u^3+lm; u-1];
        Kci = inv(Kc);
        dd = inv(Kc)*Fc*-1;
        
        
    else
        Kc = [f, 0; 0, 1];
        Fc = [u^3 ; 0];
        
        dd = inv(Kc)*Fc;
    end
    
    u = u + dd(1)
    lm = lm + dd(2)
    
    u;
end