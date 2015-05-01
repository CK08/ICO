function Data = ICOfileReader(FileName)
%
% This function reads ICO input files and creates different variables for
% the analysis cases
%
% In development
%
% Copyright (C) Diogo Cardoso, 2015
%


fid= fopen(FileName,'r');  % Open text file
line=' ';
while ~strcmp(line,'*end')
    line=fgetl(fid);
    try %Case MULTI PATCH
        if strcmp(line(1:9),'*begin_MP')
            Data = ICOreadMP(fid);
        end
    catch % Case SINGLE PATCH
        try
            if strcmp(line(1:6),'*begin')
                Data = ICOreadSP(fid,line);
            end
        catch error
            print(error);
            continue
        end
    end
end
fclose('all');




end

