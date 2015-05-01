clear all
close all
clc

scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/6.0 scrsz(4)/3.4 scrsz(3)/1.4 scrsz(4)/2.0]);

%Read KLM data --------------------------------------------------------
%File ID
fid = fopen('KLM.txt');

%Read number of patches
L = fgetl(fid);
input = textscan(L,'%d');
dim = input{1,1};

KLM = zeros(dim);

for i=1:dim
    for j=1:dim
        L = fgetl(fid);
        input = textscan(L,'%f');
        val = input{1,1};
        if(abs(val) > 1.0e-12)
            KLM(i,j) = val;
        end
    end
end

spy(KLM)

%Read Kf data ---------------------------------------------------------
%File ID
fid = fopen('Kf.txt');

%Read number of patches
L = fgetl(fid);
input = textscan(L,'%d');
dim = input{1,1};

Kf = zeros(dim);

for i=1:dim
    for j=1:dim
        L = fgetl(fid);
        input = textscan(L,'%f');
        val = input{1,1};
        if(abs(val) > 1.0e-12)
            Kf(i,j) = val;
        end
    end
end

KT = KLM;
KT(1:dim,1:dim) = KT(1:dim,1:dim) + Kf(1:dim,1:dim);

subplot(2,3,4)
spy(KLM)
grid on
%grid minor
title('KLM')

subplot(2,3,5)
spy(Kf)
grid on
%grid minor
title('Kf')

subplot(2,3,6)
spy(KT)
grid on
%grid minor
title('KLM + Kf')


%Read Kdelta data ---------------------------------------------------------
%File ID
fid = fopen('kdelta.txt');

% N = zeros(12,1);
% N(1,1) = 0.25;
% N(3,1) = 0.5;
% N(5,1) = 0.25;
% N(7,1) = -0.25;
% N(9,1) = -0.5;
% N(11,1) = -0.25;
% 
% T = zeros(12,1);
% T(2,1) = 0.25;
% T(4,1) = 0.5;
% T(6,1) = 0.25;
% T(8,1) = -0.25;
% T(10,1) = -0.5;
% T(12,1) = -0.25;
% 
% NN = T*N' + N*N' + N*T' + T*T'

%Read number of patches
KD = zeros(12);
for k1=1:3
    L = fgetl(fid);
    input = textscan(L,'%d');
    dim = input{1,1};
    
    Kdelta = zeros(dim);
    
    for i=1:dim
        for j=1:dim
            L = fgetl(fid);
            input = textscan(L,'%f');
            val = input{1,1};
            if(abs(val) > 1.0e-12)
                Kdelta(i,j) = val;
            end
        end
    end
    
    KD = KD + Kdelta;
    
    subplot(2,3,k1)
    spy(Kdelta)
    grid on
end


