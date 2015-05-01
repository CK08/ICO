function GPStress = readGPStress_MP(fid,nele1,nele2,nGP)
%
% Reads the Gauss Point Stress (GPStress.txt) file and outputs a cell for
% each patch
%
%

% Initiate variable
GP1 = [];
for i = 1:(nele1*nGP)
    GP1 = [GP1;str2num(fgetl(fid))];
end
% Alloc to structure
GPStress{1} = GP1;

% Initiate variable
GP2 = [];
for i = 1:(nele2*nGP)
    GP2 = [GP2;str2num(fgetl(fid))];
end
% Alloc to structure
GPStress{2} = GP2;


end