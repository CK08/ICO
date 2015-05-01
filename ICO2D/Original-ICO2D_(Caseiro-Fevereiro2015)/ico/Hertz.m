clear all
clc

n = 250;

E = 10.0^6;
nu = 0.3;

R = 4.0;
%x = 0.0;
%q = 0.625;

F = 96467.6; %2.0*R*q
P = F/(2.0*R);

P = F;

num = P*R*(1.0-nu^2);
den = E*pi();

b = 2.0*sqrt(num/den);
%P = 2*F/(pi()*b^2)*sqrt(b^2-x^2);

x=linspace(0,b,n);

for i=1:n
    SC(i) = 2*P/(pi()*b^2)*sqrt(b^2-x(i)^2);
end

plot(x,SC)

% clear all
% clc
% 
% syms F
% 
% E = 200.0;
% nu = 0.3;
% 
% R = 8.0;
% x = 0.0;
% 
% solve('-(E*(- x^2 - (4*F*R*(nu^2 - 1))/(E*pi))^(1/2))/(2*R*(nu^2 - 1))=9.3514', F)

% clear all
% clc
% 
% syms Sc F pi b x nu E R
% 
% num = F*R*(1.0-nu^2);
% den = E*pi;
% 
% b = 2.0*sqrt(num/den);
% 
% P = 2*F/(pi*b^2)*sqrt(b^2-x^2);
