clear all
close all
clc

%Screen data
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/2.0 scrsz(4)/2.4 scrsz(3)/2.0 scrsz(4)/2.0]);
%figure('Position',[scrsz(3)/2.03 scrsz(4)/15.0 scrsz(3)/2.0 scrsz(4)/1.175]);

%Number of divisions for each patch
np = 100;

%Activate Ligth option
Filling = 0;
Light = 1;

%View Control Net
ControlNet = 0;

%Mesh Line Options
LineWidth = 1.0;

%Contact Collocation Points
CtPoints = 0;
CPArea = 25;
RGB = [0 0 0];

%Gauss Point Data
Gauss = 1;
Field='Sxy';

%Axis
ax = [-0.1 2.0 0.0 1.1];
ax = [-0.1 1.1 -0.1 1.1];
%ax = [-0.5 12.5 0.0 6.0];
%ax = [0.0 0.3 0.0 0.3];
%ax = [-0.1 50.0 -0.1 70.0];
%ax = [-0.5 10.0 -0.5 10.0];

%Subplot Options
spdim = 2;
umesh = 1;
dmesh = 2;

%Read Initial data --------------------------------------------------------
%File ID
fid = fopen('NURBSInput.txt');

%Read number of patches
L = fgetl(fid);
input = textscan(L,'%d');
nptch = input{1,1};

%Read number of control points
ncpx = zeros(nptch,1);
ncpy = zeros(nptch,1);
for i=1:nptch
    L=fgetl(fid);
    input = textscan(L,'%d, %d');
    ncpx(i) = input{1,1};
    ncpy(i) = input{1,2};
end

%Read order of basis for each patch
p = zeros(nptch,1);
q = zeros(nptch,1);
for i=1:nptch
    L=fgetl(fid);
    input = textscan(L,'%d, %d');
    p(i) = input{1,1};
    q(i) = input{1,2};
end

%Set knot vector length
lengU = zeros(nptch,1);
lengV = zeros(nptch,1);
for i=1:nptch
    lengU(i) = ncpx(i)+p(i)+1;
    lengV(i) = ncpy(i)+q(i)+1;
end

inc = 1;
%Read control points and plot each patch
for k=1:nptch
    
    Points = zeros(ncpx(k),ncpy(k),3);
    
    %Get U knot vector
    L=fgetl(fid);
    U = str2num(L);
    
    %Get V knot vector
    L=fgetl(fid);
    V = str2num(L);
    
    %Read control points
    for i=1:ncpx(k)
        for j=1:ncpy(k)
            L=fgetl(fid);
            input = textscan(L,'%f, %f, %f');
            Points(i,j,1) = input{1,1};
            Points(i,j,2) = input{1,2};
            Points(i,j,3) = input{1,3};
%             
%             PT(inc,1) = Points(i,j,1);
%             PT(inc,2) = Points(i,j,2);
%             PT(inc,3) = Points(i,j,3);
%             
%             inc = inc + 1;
        end
    end
    
     
    for j=1:ncpy(k)
        for i=1:ncpx(k)
            
            PT(inc,1) = Points(i,j,1);
            PT(inc,2) = Points(i,j,2);
            PT(inc,3) = Points(i,j,3);
            
            inc = inc + 1;
        end
    end
    
    
    
    %Weight control points
    for k1=1:ncpx(k)
        for k2 = 1:ncpy(k)
            Pw(k1,k2,1:2) = Points(k1,k2,1:2)*Points(k1,k2,3);
            Pw(k1,k2,3)   = Points(k1,k2,3);
        end
    end
    
    %Organise data for NURBS plot
    pnts=zeros(4,ncpy(k),ncpx(k));
    for i=1:ncpx(k)
        for j=1:ncpy(k)
            pnts(1,j,i) = Pw(i,j,1);
            pnts(2,j,i) = Pw(i,j,2);
            pnts(4,j,i) = Pw(i,j,3);
        end
    end
    
    knots{1} = V;
    knots{2} = U;
    
    % Plot using nrbplot
    
    if Filling==1
        subplot(spdim,1,umesh)
        %Plot NURBS surface
        srf = nrbmak(pnts,knots);
        if(Light==1)
            nrbplot(srf,[np np],'light','on');
        else
            nrbplot(srf,[np np]);
        end
        title('Initial');
        axis auto;
        view([0 90]); % X-Y
        hold on;
        xlabel('x')
        ylabel('y')
        zlabel('z')
    end
    
    nb = lengU(k)-p-1;
    nq = lengV(k)-q-1;
    
    u = linspace(U(1),U(lengU(k)),np);
    v = linspace(V(1),V(lengV(k)),np);
    
    SRF = zeros(1,np,2);
    
    for i=p+1:lengU(k)-p
        
        k1=1;
        SRF = zeros(1,np,2);
        
        %Find span index (xi)(starting from 1)
        uspan = FindSpan(nb(k),p(k),U(i),U);
        
        %Find basis functions values (xi)
        BFu = BasisFuns(uspan,U(i),p(k),U);
        SBFu(:,k1) = BFu(:,:);
        
        for k2=1:np
            %Find span index (eta) (starting from 1)
            vspan = FindSpan(nq(k),q(k),v(k2),V);
            
            %Find basis functions values (eta)
            BFv = BasisFuns(vspan,v(k2),q(k),V);
            SBFv(:,k2) = BFv(:,:);
            
            den = 0.0;
            for l=0:q
                temp = zeros(1,1,3);
                for m=0:p
                    temp  = temp  + BFu(m+1,1)*Pw(uspan-p(k)+m, vspan-q(k)+l, :);
                    continue
                end
                
                den = den + temp(1,1,3)*BFv(l+1,1);
                SRF(k1,k2,1) = SRF(k1,k2,1) + temp(1,1,1)*BFv(l+1,1);
                SRF(k1,k2,2) = SRF(k1,k2,2) + temp(1,1,2)*BFv(l+1,1);
            end
            
            SRF(k1,k2,1) = SRF(k1,k2,1)/den;
            SRF(k1,k2,2) = SRF(k1,k2,2)/den;
            
            
        end
        
        subplot(spdim,1,umesh);
        plot(SRF(1,:,1),SRF(1,:,2),'Color','k','LineStyle','-','LineWidth',LineWidth), hold on
        title('Initial Mesh');
        axis equal;
        view([0 90]); % X-Y
        hold on;
        xlabel('x')
        ylabel('y')
        zlabel('z')
        
    end
    
    for i=q+1:lengV(k)-q
        
        k2=1;
        SRF = zeros(np,1,2);
        for k1=1:np
            
        
            %Find span index (xi)(starting from 1)
            uspan = FindSpan(nb(k),p(k),u(k1),U);

            %Find basis functions values (xi)
            BFu = BasisFuns(uspan,u(k1),p(k),U);
            SBFu(:,k1) = BFu(:,:);
        
            
            %for k2=1:np
            %Find span index (eta) (starting from 1)
            vspan = FindSpan(nq(k),q(k),V(i),V);
            
            %Find basis functions values (eta)
            BFv = BasisFuns(vspan,V(i),q(k),V);
            SBFv(:,k2) = BFv(:,:);
            
            den = 0.0;
            for l=0:q
                temp = zeros(1,1,3);
                for m=0:p
                    temp  = temp  + BFu(m+1,1)*Pw(uspan-p(k)+m, vspan-q(k)+l, :);
                    continue
                end
                
                den = den + temp(1,1,3)*BFv(l+1,1);
                SRF(k1,k2,1) = SRF(k1,k2,1) + temp(1,1,1)*BFv(l+1,1);
                SRF(k1,k2,2) = SRF(k1,k2,2) + temp(1,1,2)*BFv(l+1,1);
            end
            
            SRF(k1,k2,1) = SRF(k1,k2,1)/den;
            SRF(k1,k2,2) = SRF(k1,k2,2)/den;
            
            
        end
        
        subplot(spdim,1,umesh);
        plot(SRF(:,1,1),SRF(:,1,2),'Color','k','LineStyle','-','LineWidth',LineWidth), hold on
        title('Initial Mesh');
        axis equal;
        view([0 90]); % X-Y
        hold on;
        xlabel('x')
        ylabel('y')
        zlabel('z')
        
    end
    
    
    if(ControlNet==1)
        for k1=1:ncpx(k)
            for k2=1:ncpy(k)
                Pw(k1,k2,1:2) = Pw(k1,k2,1:2)/Pw(k1,k2,3);
            end
        end
        
        for k1=1:ncpx(k)
            for k2=1:ncpy(k)
                subplot(spdim,1,umesh)
                plot(Pw(k1,k2,1),Pw(k1,k2,2),'ko','LineWidth',1.0,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5), hold on
                axis equal;
                
                if(k1 > 1)
                    point1 = [Pw(k1-1,k2,1),Pw(k1,k2,1)];
                    point2 = [Pw(k1-1,k2,2),Pw(k1,k2,2)];
                    line(point1,point2,'Color','k','LineStyle','--','LineWidth',LineWidth), hold on
                end
                
                if(k2 > 1)
                    point1 = [Pw(k1,k2-1,1),Pw(k1,k2,1)];
                    point2 = [Pw(k1,k2-1,2),Pw(k1,k2,2)];
                    line(point1,point2,'Color','k','LineStyle','--','LineWidth',LineWidth), hold on
                end
            end
        end
    end
    
    
end

axis (ax);
%axis ([-0.5 1.5 0.0 2.0]);

fclose(fid);
%--------------------------------------------------------------------------










%Read Output data ---------------------------------------------------------
%File ID
fid = fopen('NURBSOutput.txt');

%Read number of patches
L = fgetl(fid);
input = textscan(L,'%d');
nptch = input{1,1};

%Read number of control points
ncpx = zeros(nptch,1);
ncpy = zeros(nptch,1);
for i=1:nptch
    L=fgetl(fid);
    input = textscan(L,'%d, %d');
    ncpx(i) = input{1,1};
    ncpy(i) = input{1,2};
end

%Read order of basis for each patch
p = zeros(nptch,1);
q = zeros(nptch,1);
for i=1:nptch
    L=fgetl(fid);
    input = textscan(L,'%d, %d');
    p(i) = input{1,1};
    q(i) = input{1,2};
end

%Set knot vector length
lengU = zeros(nptch,1);
lengV = zeros(nptch,1);
for i=1:nptch
    lengU(i) = ncpx(i)+p(i)+1;
    lengV(i) = ncpy(i)+q(i)+1;
end

%Read control points and plot each patch
inc = 1;
for k=1:nptch
    
    Points = zeros(ncpx(k),ncpy(k),3);
    
    %Get U knot vector
    L=fgetl(fid);
    U = str2num(L);
    
    %Get V knot vector
    L=fgetl(fid);
    V = str2num(L);
    
    %Read control points
    for i=1:ncpx(k)
        for j=1:ncpy(k)
            L=fgetl(fid);
            input = textscan(L,'%f, %f, %f');
            Points(i,j,1) = input{1,1};
            Points(i,j,2) = input{1,2};
            Points(i,j,3) = input{1,3};
        end
    end
    
    for j=1:ncpy(k)
        for i=1:ncpx(k)
            
            PT_def(inc,1) = Points(i,j,1);
            PT_def(inc,2) = Points(i,j,2);
            PT_def(inc,3) = Points(i,j,3);
            
            inc = inc + 1;
        end
    end
    
    %Weight control points
    for k1=1:ncpx(k)
        for k2 = 1:ncpy(k)
            Pw(k1,k2,1:2) = Points(k1,k2,1:2)*Points(k1,k2,3);
            Pw(k1,k2,3)   = Points(k1,k2,3);
        end
    end
    
    %Organise data for NURBS plor
    pnts=zeros(4,ncpy(k),ncpx(k));
    for i=1:ncpx(k)
        for j=1:ncpy(k)
            pnts(1,j,i) = Pw(i,j,1);
            pnts(2,j,i) = Pw(i,j,2);
            pnts(4,j,i) = Pw(i,j,3);
        end
    end
    
    knots{1} = V;
    knots{2} = U;
    
    %  Plot using nrbplot    
    
    if Filling==1
        subplot(spdim,1,dmesh)
        %Plot NURBS surface
        srf = nrbmak(pnts,knots);
        if(Light==1)
            nrbplot(srf,[np np],'light','on');
        else
            nrbplot(srf,[np np]);
        end
        title('Deformed Mesh');
        axis equal;
        view([0 90]); % X-Y
        hold on;
        xlabel('x')
        ylabel('y')
        zlabel('z')
    end
    
    nb = lengU(k)-p-1;
    nq = lengV(k)-q-1;
    
    u = linspace(U(1),U(lengU(k)),np);
    v = linspace(V(1),V(lengV(k)),np);
    
    SRF = zeros(1,np,2);
    
    for i=p+1:lengU(k)-p
        
        k1=1;
        SRF = zeros(1,np,2);
        
        %Find span index (xi)(starting from 1)
        uspan = FindSpan(nb(k),p(k),U(i),U);
        
        %Find basis functions values (xi)
        BFu = BasisFuns(uspan,U(i),p(k),U);
        SBFu(:,k1) = BFu(:,:);
        
        for k2=1:np
            %Find span index (eta) (starting from 1)
            vspan = FindSpan(nq(k),q(k),v(k2),V);
            
            %Find basis functions values (eta)
            BFv = BasisFuns(vspan,v(k2),q(k),V);
            SBFv(:,k2) = BFv(:,:);
            
            den = 0.0;
            for l=0:q
                temp = zeros(1,1,3);
                for m=0:p
                    temp  = temp  + BFu(m+1,1)*Pw(uspan-p(k)+m, vspan-q(k)+l, :);
                    continue
                end
                
                den = den + temp(1,1,3)*BFv(l+1,1);
                SRF(k1,k2,1) = SRF(k1,k2,1) + temp(1,1,1)*BFv(l+1,1);
                SRF(k1,k2,2) = SRF(k1,k2,2) + temp(1,1,2)*BFv(l+1,1);
            end
            
            SRF(k1,k2,1) = SRF(k1,k2,1)/den;
            SRF(k1,k2,2) = SRF(k1,k2,2)/den;
            
            
        end
        
        subplot(spdim,1,dmesh);
        plot(SRF(1,:,1),SRF(1,:,2),'Color','k','LineStyle','-','LineWidth',LineWidth), hold on
        title('Deformed Mesh');
        axis equal;
        view([0 90]); % X-Y
        hold on;
        xlabel('x')
        ylabel('y')
        zlabel('z')
        
    end
    
    
    for i=q+1:lengV(k)-q
        
        k2=1;
        SRF = zeros(np,1,2);
        for k1=1:np
            
        
            %Find span index (xi)(starting from 1)
            uspan = FindSpan(nb(k),p(k),u(k1),U);

            %Find basis functions values (xi)
            BFu = BasisFuns(uspan,u(k1),p(k),U);
            SBFu(:,k1) = BFu(:,:);
        
            
            %for k2=1:np
            %Find span index (eta) (starting from 1)
            vspan = FindSpan(nq(k),q(k),V(i),V);
            
            %Find basis functions values (eta)
            BFv = BasisFuns(vspan,V(i),q(k),V);
            SBFv(:,k2) = BFv(:,:);
            
            den = 0.0;
            for l=0:q
                temp = zeros(1,1,3);
                for m=0:p
                    temp  = temp  + BFu(m+1,1)*Pw(uspan-p(k)+m, vspan-q(k)+l, :);
                    continue
                end
                
                den = den + temp(1,1,3)*BFv(l+1,1);
                SRF(k1,k2,1) = SRF(k1,k2,1) + temp(1,1,1)*BFv(l+1,1);
                SRF(k1,k2,2) = SRF(k1,k2,2) + temp(1,1,2)*BFv(l+1,1);
            end
            
            SRF(k1,k2,1) = SRF(k1,k2,1)/den;
            SRF(k1,k2,2) = SRF(k1,k2,2)/den;
            
            
        end
        
        subplot(spdim,1,dmesh);
        plot(SRF(:,1,1),SRF(:,1,2),'Color','k','LineStyle','-','LineWidth',LineWidth), hold on
        title('Deformed Mesh');
        axis equal;
        view([0 90]); % X-Y
        hold on;
        xlabel('x')
        ylabel('y')
        zlabel('z')
        
    end
    
    if(ControlNet==1)
        for k1=1:ncpx(k)
            for k2=1:ncpy(k)
                Pw(k1,k2,1:2) = Pw(k1,k2,1:2)/Pw(k1,k2,3);
            end
        end
        
        for k1=1:ncpx(k)
            for k2=1:ncpy(k)
                subplot(spdim,1,dmesh)
                plot(Pw(k1,k2,1),Pw(k1,k2,2),'ko','LineWidth',1.0,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5), hold on
                axis equal;
                view([0 90]); % X-Y
                hold on;
                xlabel('x')
                ylabel('y')
                zlabel('z')
                
                if(k1 > 1)
                    point1 = [Pw(k1-1,k2,1),Pw(k1,k2,1)];
                    point2 = [Pw(k1-1,k2,2),Pw(k1,k2,2)];
                    line(point1,point2,'Color','k','LineStyle','--','LineWidth',LineWidth), hold on
                end
                
                if(k2 > 1)
                    point1 = [Pw(k1,k2-1,1),Pw(k1,k2,1)];
                    point2 = [Pw(k1,k2-1,2),Pw(k1,k2,2)];
                    line(point1,point2,'Color','k','LineStyle','--','LineWidth',LineWidth), hold on
                end
            end
        end
    end
    
end

axis (ax);

if(CtPoints == 1)
    
    %Read contact data
    fid = fopen('ContactData.txt');
    L = fgetl(fid);
    input = textscan(L,'%d');
    p_s = input{1,1};
    
    L = fgetl(fid);
    input = textscan(L,'%d');
    ncp_s = input{1,1};
    
    L=fgetl(fid);
    U_s = str2num(L);
    
    L=fgetl(fid);
    conn_s = str2num(L);
    
    %Compute greville points
    grev=zeros(ncp_s,1);
    for i=1:ncp_s
        for j=1:p_s
            grev(i) = grev(i) + U_s(i+j)/(double(p_s));
        end
    end
    
    %Undeformed configuration
    for k1=1:length(conn_s)
        Pt_sl(k1,1:3) = PT(conn_s(k1),1:3);
    end
    
        
    nbs = length(U_s)-p_s-1;
    Curve = zeros(ncp_s,2);
    
    for k1=1:ncp_s
        %Find span index (starting from 1)
        s = FindSpan(nbs,p_s,grev(k1),U_s);
        
        %Find basis functions values
        BF = BasisFuns(s,grev(k1),p_s,U_s);
        
        Y(:,k1) = BF(:,:);
        
        den = 0.0;
        for i=0:p_s
            den = den + Y(i+1,k1)*Pt_sl(s-p_s+i,3);
            Curve(k1,1:2) = Curve(k1,1:2) + Y(i+1,k1)*Pt_sl(s-p_s+i,1:2)*Pt_sl(s-p_s+i,3);
        end
        
        Curve(k1,:) = Curve(k1,:)/den;
        
%         if (k1 == 1)
%             scur = s;
%         else
%             if(s == scur +1)
%                 scur = s;
%                 Node(nnode,1) = Curve(k1,1);
%                 Node(nnode,2) = Curve(k1,2);
%                 nnode = nnode + 1;
%             end
%         end
        
        subplot(spdim,1,umesh)
        scatter(Curve(k1,1),Curve(k1,2),CPArea,RGB,'fill','d'), hold on

    end
    
    %Deformed configuration
    for k1=1:length(conn_s)
        Pt_sl(k1,1:3) = PT_def(conn_s(k1),1:3);
    end
    
        
    nbs = length(U_s)-p_s-1;
    Curve = zeros(ncp_s,2);
    
    for k1=1:ncp_s
        %Find span index (starting from 1)
        s = FindSpan(nbs,p_s,grev(k1),U_s);
        
        %Find basis functions values
        BF = BasisFuns(s,grev(k1),p_s,U_s);
        
        Y(:,k1) = BF(:,:);
        
        den = 0.0;
        for i=0:p_s
            den = den + Y(i+1,k1)*Pt_sl(s-p_s+i,3);
            Curve(k1,1:2) = Curve(k1,1:2) + Y(i+1,k1)*Pt_sl(s-p_s+i,1:2)*Pt_sl(s-p_s+i,3);
        end
        
        Curve(k1,:) = Curve(k1,:)/den;

        subplot(spdim,1,dmesh)
        scatter(Curve(k1,1),Curve(k1,2),CPArea,RGB,'fill','d'), hold on

    end
 
end

fclose(fid);
%--------------------------------------------------------------------------


if(Gauss == 1)
    fid = fopen('GPCoords.txt');
    L = fgetl(fid);
    input = textscan(L,'%d');
    icycle = input{1,1};
    
    xGP = zeros(icycle,1);
    yGP = zeros(icycle,1);
    
    for i=1:icycle
        L = fgetl(fid);
        input = textscan(L,'%f, %f');
        xGP(i) = input{1,1};
        yGP(i) = input{1,2};
    end

    fid = fopen('GPStress.txt');
    Sxx = zeros(icycle,1);
    Syy = zeros(icycle,1);
    Sxy = zeros(icycle,1);
    
    for i=1:icycle
        L = fgetl(fid);
        input = textscan(L,'%f, %f, %f');
        Sxx(i) = input{1,1};
        Syy(i) = input{1,2};
        Sxy(i) = input{1,3};
    end
    
    fid = fopen('GPStrain.txt');
    Exx = zeros(icycle,1);
    Eyy = zeros(icycle,1);
    Exy = zeros(icycle,1);
    
    for i=1:icycle
        L = fgetl(fid);
        input = textscan(L,'%f, %f, %f');
        Exx(i) = input{1,1};
        Eyy(i) = input{1,2};
        Exy(i) = input{1,3};
    end
    
    if(Field == 'Sxx')
        data = Sxx;
    elseif(Field == 'Syy')
        data = Syy;
    elseif(Field == 'Sxy')
        data = Sxy;
    elseif(Field == 'Exx')
        data = Exx;    
    elseif(Field == 'Eyy')
        data = Eyy;    
    elseif(Field == 'Exy')
        data = Exy;
    end    
        
        
    subplot(spdim,1,dmesh)
    scatter(xGP,yGP,15,data,'fill','s')
    t = colorbar('peer',gca);
    set(get(t,'xlabel'),'String', Field);
    %colorbar
    
%     figure
%     %plot(xGP,yGP,'.','markersize',12)
%     [xq,yq] = meshgrid(xGP,yGP);
%     vq = griddata(xGP,yGP,data,xq,yq);
%     surf(xq,yq,vq)
%     colorbar
    
end


% %Read stress in the undeformed mesh
% %File ID
% fid = fopen('StressI2.txt');
% 
% %Read number of patches
% L = fgetl(fid);
% input = textscan(L,'%d');
% nstr = input{1,1};
% 
% for i=1:nstr
%     L=fgetl(fid);
%     data(i,:) = str2num(L);
% end
% 
% for i=1:nstr
%     x(i,1) = data(i,1);
%     y(i,1) = data(i,2);
%     s1(i,1) = data(i,3);
% end
% 
% F = TriScatteredInterp(x(:),y(:),s1(:));
% [X,Y] = meshgrid(x,y);
% S1 = F(X,Y);
% 
% % subplot(3,1,3)
% % scatter3(x,y,s1);
% % view([0 90]);
% % mesh(X,Y,S1)
% % axis auto;
% 
% subplot(3,1,1)
% %contour(X,Y,S1,'ShowText','on');
% 
% contourf(X,Y,S1,10);
% 
% axis equal;
% view([0 90]); % X-Y
% hold on;
% xlabel('x')
% ylabel('y')
% zlabel('z')
% 
% fclose(fid);

