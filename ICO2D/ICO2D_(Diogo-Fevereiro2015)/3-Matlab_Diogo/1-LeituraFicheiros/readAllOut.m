function geometry = readAllOut(PathAllOut, PathInput)
%----------------------------------------------------------------
%
% Read function for the NURBSAllOut.txt file from ICO2D
%
% Generates a structure called geometry with number fields from
% zero (initial state) to the maximum of increments or even 
% iteration saved in the NURBSAllOut.txt
%
% Copyright (C) Diogo Cardoso, 2015
%
%----------------------------------------------------------------


% Get the number of patches
fid3 = fopen(PathAllOut,'r');
npatch = str2double(fgetl(fid3));
fclose(fid3);
% Initial structure
fid2 = fopen(PathInput,'r');
[ncpi, ordersi, knotsi, coefsi] = readNURBS_MP(fid2);
geometry.('f0').coefs = coefsi;
geometry.('f0').order = ordersi;
geometry.('f0').number = ncpi;
geometry.('f0').knots = knotsi;
fclose(fid2);
% Analysis structure
fid = fopen(PathAllOut,'r');
pos = 0;
if npatch>=2 %MultiplePatches
    while 1
        try
            [ncpf, ordersf, knotsf, coefsf] = readNURBS_MP(fid);
            pos = pos + 1; 
            geometry.(['f' num2str(pos)]).coefs = coefsf;
            geometry.(['f' num2str(pos)]).order = ordersf;
            geometry.(['f' num2str(pos)]).number = ncpf;
            geometry.(['f' num2str(pos)]).knots = knotsf;
        catch
            break;
        end
    end
    fclose(fid);
    
else %SinglePatch
    %ONGOING....................
    geometry = 0;
end

end