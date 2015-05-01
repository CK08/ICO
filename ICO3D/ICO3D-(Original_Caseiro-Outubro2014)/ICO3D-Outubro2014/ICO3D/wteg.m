clear all
clc

B = bucky;
r = symrcm(B);
m = symamd(B);

figure
subplot(1,3,1)
spy(B)
title('Original')

subplot(1,3,2)
spy(B(r,r))
title('Reverse Cuthill-McKee')

subplot(1,3,3)
spy(B(m,m))
title('Approx Min Degree')