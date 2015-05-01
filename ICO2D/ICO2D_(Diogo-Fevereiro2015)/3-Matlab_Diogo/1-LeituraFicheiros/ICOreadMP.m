function Data = ICOreadMP(fid)
%
%
%
line=' ';
while ~strcmp(line,'*end')
    line=fgetl(fid);
    %----------------------------------------------------------------------
    % Case *begins
    %----------------------------------------------------------------------
    try
        if strcmp(line(1:6),'*begin')
            npatch=fgetl(fid);npatch=str2double(npatch(1)); % npatch of the problem
            Data.npatch = npatch;
            dim=fgetl(fid);dimension=str2double(dim(1)); % Dimension of the problem
            if dimension==3,dd=fgetl(fid);direction_degrees=str2num(dd(1:7));end % Degree in each direction
            if dimension==2,dd=fgetl(fid);direction_degrees=str2num(dd(1:4));end % Degree in each direction
            Data.patch1.order = direction_degrees;
            if dimension==3,dd=fgetl(fid);direction_degrees=str2num(dd(1:7));end % Degree in each direction
            if dimension==2,dd=fgetl(fid);direction_degrees=str2num(dd(1:4));end % Degree in each direction
            Data.patch2.order = direction_degrees;
            if dimension==3,ncp=fgetl(fid);number_cp=str2num(ncp(1:9));end % Number of Control Points
            if dimension==2,ncp=fgetl(fid);number_cp=str2num(ncp(1:7));end % Number of Control Points
            Data.patch1.number = number_cp;
            if dimension==3,ncp=fgetl(fid);number_cp=str2num(ncp(1:9));end % Number of Control Points
            if dimension==2,ncp=fgetl(fid);number_cp=str2num(ncp(1:7));end % Number of Control Points
            Data.patch2.number = number_cp;
            if dimension==3,tkv=fgetl(fid);knot_vector_type=str2num(tkv(1:7));end % Knot Vector Type
            if dimension==2,tkv=fgetl(fid);knot_vector_type=str2num(tkv(1:4));end % Knot Vector Type
            if dimension==3,tkv=fgetl(fid);knot_vector_type=str2num(tkv(1:7));end % Knot Vector Type
            if dimension==2,tkv=fgetl(fid);knot_vector_type=str2num(tkv(1:4));end % Knot Vector Type
            ncp_total = fgetl(fid);Data.ncp_total = str2num(ncp_total(1:3));
            clear dim dd ncp tkv;
        end
    catch
    end
    %----------------------------------------------------------------------
    % Case knots
    %----------------------------------------------------------------------
    try
        
        if strcmp(line(1:6),'*knots')
            if dimension==3                     % Knot vectors for each direction
                kv11=fgetl(fid);KV11=str2num(kv11);
                kv21=fgetl(fid);KV21=str2num(kv21);
                kv31=fgetl(fid);KV31=str2num(kv31);
                kv12=fgetl(fid);KV12=str2num(kv12);
                kv22=fgetl(fid);KV22=str2num(kv22);
                kv32=fgetl(fid);KV32=str2num(kv32);
                Data.patch1.knots = {KV11,KV21,KV31};
                Data.patch2.knots = {KV12,KV22,KV32};
                clear kv11 kv21 kv31 kv12 kv22 kv32;
            elseif dimension==2
                kv11=fgetl(fid);KV11=str2num(kv11);
                kv21=fgetl(fid);KV21=str2num(kv21);
                kv12=fgetl(fid);KV12=str2num(kv12);
                kv22=fgetl(fid);KV22=str2num(kv22);
                Data.patch1.knots = {KV11,KV21};
                Data.patch2.knots = {KV12,KV22};
                clear kv11 kv21 kv12 kv22;
            end
        end
    catch
    end
    %----------------------------------------------------------------------
    % Case bnet
    %----------------------------------------------------------------------
    try
        if strcmp(line(1:5),'*bnet')
            if dimension==3
                number1 = Data.patch1.number;
                coefs1 = zeros(4,number1(1),number1(2),number1(3));
                for i = 1:number1(1)
                    for j = 1:number1(2)
                        for k = 1:number1(3)
                            coor = str2num(fgetl(fid));
                            coefs1(1,i,j,k) = coor(1)*coor(4);
                            coefs1(2,i,j,k) = coor(2)*coor(4);
                            coefs1(3,i,j,k) = coor(3)*coor(4);
                            coefs1(4,i,j,k) = coor(4);
                        end
                    end
                end
                
                number2 = Data.patch2.number;
                coefs2 = zeros(4,number2(1),number2(2),number2(3));
                
                for i = 1:number2(1)
                    for j = 1:number2(2)
                        for k = 1:number2(3)
                            coor = str2num(fgetl(fid));
                            coefs2(1,i,j,k) = coor(1)*coor(4);
                            coefs2(2,i,j,k) = coor(2)*coor(4);
                            coefs2(3,i,j,k) = coor(3)*coor(4);
                            coefs2(4,i,j,k) = coor(4);
                        end
                    end
                end
                Data.patch1.coefs = coefs1;
                Data.patch2.coefs = coefs2;
            end
            if dimension==2
                number1 = Data.patch1.number;
                coefs1 = zeros(3,number1(1),number1(2));
                for i = 1:number1(1)
                    for j = 1:number1(2)
                        coor = str2num(fgetl(fid));
                        coefs1(1,i,j) = coor(1)*coor(3);
                        coefs1(2,i,j) = coor(2)*coor(3);
                        coefs1(3,i,j) = coor(3);
                    end
                end
                
                number2 = Data.patch2.number;
                coefs2 = zeros(3,number2(1),number2(2));
                
                for i = 1:number2(1)
                    for j = 1:number2(2)
                        coor = str2num(fgetl(fid));
                        coefs2(1,i,j) = coor(1)*coor(3);
                        coefs2(2,i,j) = coor(2)*coor(3);
                        coefs2(3,i,j) = coor(3);
                    end
                end
                Data.patch1.coefs = coefs1;
                Data.patch2.coefs = coefs2;
                clear number1 number2 coefs1 coefs2 i j coor;
            end
        end
    catch
    end
    %----------------------------------------------------------------------
    % case *Multistep
    %----------------------------------------------------------------------
    try
        if strcmp(line(1:9),'*Multistep')
            nsteps = fgetl(fid); nsteps = str2double(nsteps);
            Data.Multistep.nsteps = nsteps;
            for step_cycle=1:nsteps
                line = fgetl(fid);
                %----------------------------------------------------------
                % case *Step
                %----------------------------------------------------------
                if strcmp(line(1:5),'*step')
                    while ~strcmp(line,'*endstep')
                    line = fgetl(fid);
                    %------------------------------------------------------
                    % case *Increment
                    %------------------------------------------------------
                    try
                        if strcmp(line(1:10),'*increment')
                            increment = fgetl(fid);
                            increment = str2double(increment);
                            Data.Multistep.(['step' num2str(step_cycle)]).increment =...
                                increment;
                        end
                    catch
                    end
                    %------------------------------------------------------
                    % case *Iteration
                    %------------------------------------------------------
                    try
                        if strcmp(line(1:10),'*iteration')
                            iteration = fgetl(fid);
                            increiterationment = str2double(iteration);
                            Data.Multistep.(['step' num2str(step_cycle)]).iteration =...
                                iteration;
                        end
                    catch
                    end
                    %------------------------------------------------------
                    % case *BCdof
                    %------------------------------------------------------
                    
                    %------------------------------------------------------
                    % case *DISPdof
                    %------------------------------------------------------
                    
                    %------------------------------------------------------
                    % case *Iteration
                    %------------------------------------------------------
                    end
                end
            end
            clear step_cycle    
        end
        
    catch
    end
    %---------------------------------------------------------------------- 
    % MP_connectivity
    %----------------------------------------------------------------------
    try
        if strcmp(line(1:8),'*MP_conn')
            nele = (sum(ceil(diff(Data.patch1.knots{1})))*sum(ceil(diff(Data.patch1.knots{2}))))...
                + (sum(ceil(diff(Data.patch2.knots{1})))*sum(ceil(diff(Data.patch2.knots{2}))));
            conn = [];
            for i = 1:nele
                ele = fgetl(fid); ele = str2num(ele);
                conn(i,:)=ele;
            end
            clear ele nele;
            Data.conn = conn;
        end
    catch
    end
end

end