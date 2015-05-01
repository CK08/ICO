function Data = ICOreadSP(fid,line)
%
%
%

while ~strcmp(line,'*end')
    
    try
        % Case *begins
        if strcmp(line(1:6),'*begin')
            dim=fgetl(fid);dimension=str2double(dim(1)); % Dimension of the problem
            if dimension==3,dd=fgetl(fid);direction_degrees=str2num(dd(1:7));end % Degree in each direction
            if dimension==2,dd=fgetl(fid);direction_degrees=str2num(dd(1:4));end % Degree in each direction
            Data.patch1.order = direction_degrees;
            if dimension==3,ncp=fgetl(fid);number_cp=str2num(ncp(1:7));end % Number of Control Points
            if dimension==2,ncp=fgetl(fid);number_cp=str2num(ncp(1:4));end % Number of Control Points
            Data.patch1.number = number_cp;
            if dimension==3,tkv=fgetl(fid);knot_vector_type=str2num(tkv(1:7));end % Knot Vector Type
            if dimension==2,tkv=fgetl(fid);knot_vector_type=str2num(tkv(1:4));end % Knot Vector Type
            clear dim dd ncp tkv;
        end
        
        % Case knots
        if strcmp(line(1:6),'*knots')
            if dimension==3                     % Knot vectors for each direction
                kv1=fgetl(fid);KV1=str2num(kv1);Data.patch1.knots{1} = kv1;
                kv2=fgetl(fid);KV2=str2num(kv2);Data.patch1.knots{2} = kv2;
                kv3=fgetl(fid);KV3=str2num(kv3);Data.patch1.knots{3} = kv3;
                clear kv1 kv2 kv3;
            end
            if dimension==2
                kv1=fgetl(fid);KV1=str2num(kv1);Data.patch1.knots{1} = kv1;
                kv2=fgetl(fid);KV2=str2num(kv2);Data.patch1.knots{2} = kv2;
                clear kv1 kv2;
            end
        end
        % Case Control Points
        if strcmp(line(1:5),'*bnet')
            % coefs initialization
            coefs = [];
            if dimension==3
                for i=1:(number_cp(1)*number_cp(2)*number_cp(3))
                    cp=fgetl(fid);
                    cp=str2num(cp);
                    coefs = [coefs;cp];
%                     for j=1:4
%                         CP{i,j}=num2str(cp(j),7);
%                     end
                end
            elseif dimension==2
                for i=1:(number_cp(1)*number_cp(2))
                    cp=fgetl(fid);
                    cp=str2num(cp);
                    coefs = [coefs;cp];
%                     for j=1:3   
%                         % OPTION FOR X,Y,Z,W coordinates for 2-D
%                         if j==3
%                             CP{i,3}='0';
%                             CP{i,4}=num2str(cp(j),7);
%                         else
%                             CP{i,j}=num2str(cp(j),7);
%                         end
%                     end
                end
            end
            Data.patch1.coefs = coefs;
            clear ans cp;
        end
        
        % case *element
        if strcmp(line(1:8),'*element')
            element = fgetl(fid);
            Data.patch1.element = element;
            clear element;
        end
        
        % case *material
        if strcmp(line(1:9),'*material')
            nprop = fgetl(fid); nprop = str2double(nprop);
            prop = fgetl(fid); prop = str2num(prop);
            Data.patch1.material.nprop = nprop;
            Data.patch1.material.prop = prop;
            clear nprop prop;
        end
        
        % case *Multistep
        if strcmp(line(1:9),'*Multistep')
            nsteps = fgetl(fid); nsteps = str2double(nsteps);
            Data.Multistep.nsteps = nsteps;
            
            
        end
        
        
        
        
        
        
        
        
        
        line=fgetl(fid);
    catch
        line=fgetl(fid);
        continue;
    end
end

end