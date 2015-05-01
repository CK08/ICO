function arrows = nrbLoadplotx(nurbs, data, str)
%
% The nrbLoadplotx function creates the representation of the nurbs
% structure load in the current axes.
%
% Input:  - nurbs structure (see nrbmak)
%         - data: [n x 3] matrix with the direction of the load
%         to apply in the second column, magnitude on the third and 
%         the corresponding index point on the first;
%         - color (RGB vector color)
%
% Output: - [inic, end] of the arrows
%
% Copyright (C) Diogo Cardoso 2014
%
if (~isstruct(nurbs)) || (size(data,2)~=3),return;end

coor = nrbind2coor(nurbs,data(:,1));

if max(size(nurbs.order))==3
    np = [100 100 100];
elseif max(size(nurbs.order))==2
    np = [100 100];
else
    return;
end

nrbplotx(nurbs,np);hold on;
arrows = LoadArrow(coor,data,str);
view([1 1 1]);
end

function arrows = LoadArrow(PointCoor,data,str)
%
%
%
pos = 1;
for i =1:size(data,1)
    if data(i,2) == (1) && data(i,3)>0
        inic = [PointCoor(i,1)-1,PointCoor(i,2),PointCoor(i,3)];
        fim  = [PointCoor(i,1),PointCoor(i,2),PointCoor(i,3)];
    elseif data(i,2) == (1) && data(i,3)<0
        inic = [PointCoor(i,1),PointCoor(i,2),PointCoor(i,3)];
        fim  = [PointCoor(i,1)-1,PointCoor(i,2),PointCoor(i,3)];
    elseif data(i,2) == (2) && data(i,3)>0
        inic = [PointCoor(i,1),PointCoor(i,2)-1,PointCoor(i,3)];
        fim  = [PointCoor(i,1),PointCoor(i,2),PointCoor(i,3)];
    elseif data(i,2) == (2) && data(i,3)<0
        inic = [PointCoor(i,1),PointCoor(i,2),PointCoor(i,3)];
        fim  = [PointCoor(i,1),PointCoor(i,2)-1,PointCoor(i,3)];
    elseif data(i,2) == (3) && data(i,3)>0
        inic = [PointCoor(i,1),PointCoor(i,2),PointCoor(i,3)-1];
        fim  = [PointCoor(i,1),PointCoor(i,2),PointCoor(i,3)];
    elseif data(i,2) == (3) && data(i,3)<0
        inic = [PointCoor(i,1),PointCoor(i,2),PointCoor(i,3)+1];
        fim  = [PointCoor(i,1),PointCoor(i,2),PointCoor(i,3)];
    else
        error('Error: unknowned direction value');
    end
    arrows(pos,:) = [inic, fim];
    pos = pos + 1;
    arrow3d(inic,fim,15,'cylinder',[0.5,0.5]);
    hold on;
    text(inic(1)+0.1,inic(2)+0.1,inic(3)+0.1,[num2str(data(i,3)) str],'FontSize',14);
end
end