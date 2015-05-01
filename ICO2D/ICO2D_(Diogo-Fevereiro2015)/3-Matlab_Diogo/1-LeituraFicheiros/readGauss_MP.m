function GP = readGauss_MP(fid,nele1,nele2)
%
%
%
%

total_GP = str2double(fgetl(fid));
nGP = total_GP/(nele1+nele2);
% Initiate variable
GP1 = [];
for i = 1:(nele1*nGP)
    GP1 = [GP1;str2num(fgetl(fid))];
end
% Alloc to structure
GP{1} = GP1;

% Initiate variable
GP2 = [];
for i = 1:(nele2*nGP)
    GP2 = [GP2;str2num(fgetl(fid))];
end
% Alloc to structure
GP{2} = GP2;


end