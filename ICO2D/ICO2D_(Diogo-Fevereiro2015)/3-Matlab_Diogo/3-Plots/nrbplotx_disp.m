function nrbplotx_disp(nurbs_1, nurbs_2, subd, dir)
%
% 
% NRBPLOTx_disp: Plot a NURBS surface, or the boundary of a NURBS volume
%           with the respective span limits and displacement contour 
%           representation.
% 
% Calling Sequence:
% 
%   nrbplotx_disp (nrbs1, nrbs2, subd, dir)
% 
% INPUT:
% 
%   nrbs_1  : NURBS surface or volume for the first increment, see nrbmak.
% 
%   nrbs_2  : NURBS surface or volume for the second increment, see nrbmak.
%
%   npnts	: Number of evaluation points, for a surface or volume, a row 
%       vector with the number of points along each direction.
% 
%   dir     : Direction of the displacement representation.
%               1 - xx direction
%               2 - yy direction
%               3 - zz direction
%               4 - overall
%
%
%   Based on nrbplot from:    
%   -----------------------------------------------------------------------
%    Copyright (C) 2000 Mark Spink
%    Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
%    Copyright (C) 2012 Rafael Vazquez
%   -----------------------------------------------------------------------
%
%
%    Copyright (C) 2014 Diogo Cardoso
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
nargs = nargin;
if nargs < 2
  error ('Not enought parameters');
elseif nargs >4
  error ('Too much parameters');
end
% teste!!!!!!!!!!!!!
colormap('default')

% convert the number of subdivisions in number of points
subd = subd+1;

% plot the curve or surface
if (iscell (nurbs_1.knots))
 if (size (nurbs_1.knots,2) == 2) % plot a NURBS surface
  knt = nurbs_1.knots;
  p = nrbeval (nurbs_1, {linspace(knt{1}(1),knt{1}(end),subd(1)) ...
                       linspace(knt{2}(1),knt{2}(end),subd(2))});
  xx = squeeze(p(1,:,:));
  yy = squeeze(p(2,:,:));
  zz = squeeze(p(3,:,:));
  
  knt2 = nurbs_2.knots; 
  p2 = nrbeval (nurbs_2, {linspace(knt2{1}(1),knt2{1}(end),subd(1)) ...
                       linspace(knt2{2}(1),knt2{2}(end),subd(2))});
  xx2 = squeeze(p2(1,:,:));
  yy2 = squeeze(p2(2,:,:));
  zz2 = squeeze(p2(3,:,:));
  
  if dir==1 % xx disp
      disp = xx-xx2;
  elseif dir==2 % yy disp
      disp = yy-yy2; 
  elseif dir==3 % zz disp
      disp = zz-zz2;
  elseif dir==4 % overall disp
      disp = sqrt((xx-xx2).^2+(yy-yy2).^2+(zz-zz2).^2);
  end
  surf(xx2,yy2,zz2,disp,'EdgeColor','none');hold on;
  %colorbar;
  plotSpanLim(nurbs_2(1),subd(2));
  
 elseif (size (nurbs_1.knots,2) == 3) % plot the boundaries of a NURBS volume
  bnd = nrbextract (nurbs_1);
  bnd2 = nrbextract (nurbs_2);
  hold_flag = ishold;
  hold on;
  nrbplotx_disp (bnd(1), bnd2(1), subd(2:3), dir);
  nrbplotx_disp (bnd(2), bnd2(2), subd(2:3), dir);
  nrbplotx_disp (bnd(3), bnd2(3), subd(2:3), dir);
  nrbplotx_disp (bnd(4), bnd2(4), subd(2:3), dir);
  nrbplotx_disp (bnd(5), bnd2(5), subd(2:3), dir);
  nrbplotx_disp (bnd(6), bnd2(6), subd(2:3), dir);
%   colorbar;
%   plotSpanLim(bnd2(1),subd(2));
%   plotSpanLim(bnd2(2),subd(2));
%   plotSpanLim(bnd2(3),subd(1));
%   plotSpanLim(bnd2(4),subd(1));
%   plotSpanLim(bnd2(5),subd(1));
%   plotSpanLim(bnd2(6),subd(1));

  if (~hold_flag)
    hold off
  end
 
 else
  error ('nrbplotx_disp: some argument is not correct')
 end
else
  error('Error: No nrbplotx_disp for curves');
end
%axis equal;

end