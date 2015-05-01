clear all
clc

syms dx dy Rx xi yi

dx = Rx*xi
dy = Rx*yi

f = sqrt(dx*dx+dy*dy)

simplify(collect(f,Rx))

a=[1 2 3 4 5];
b = [1.23 3.3465 5.123 0.3456 1.7456]

res1 = 0;
res2 = 0;
for i=1:5
    res1 = res1 + a(i)*a(i)*b(i)*b(i);
    res2 = res2 + (a(i)*b(i))^2;
end

res1 
res2