function nrbplotx(nurbs, subd, varargin)
% 
% NRBPLOTx: Plot a NURBS curve or surface, or the boundary of a NURBS volume
%           with the respective span limits.
% 
% Calling Sequence:
% 
%   nrbplot (nrb, subd)
%   nrbplot (nrb, subd, p, v)
% 
% INPUT:
% 
%   nrb		: NURBS curve, surface or volume, see nrbmak.
% 
%   npnts	: Number of evaluation points, for a surface or volume, a row 
%       vector with the number of points along each direction.
% 
%   [p,v]       : property/value options
%
%               Valid property/value pairs include:
%
%               Property        Value/{Default}
%               -----------------------------------
%               light           {off} | on
%               colormap        {'copper'}
%
% Example:
%
%   Plot the test surface with 20 points along the U direction
%   and 30 along the V direction
%
%   nrbplot(nrbtestsrf, [20 30])
%
%    Copyright (C) 2000 Mark Spink
%    Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
%    Copyright (C) 2012 Rafael Vazquez
%    Copyright (C) 2014 Diogo Cardoso
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
nargs = nargin;
if nargs < 2
  error ('Need a NURBS to plot and the number of subdivisions!');
elseif rem(nargs+2,2)
  error ('Param value pairs expected')
end

% Default values
light='off';
cmap='summer';

% Recover Param/Value pairs from argument list
for i=1:2:nargs-2
  Param = varargin{i};
  Value = varargin{i+1};
  if (~ischar (Param))
    error ('Parameter must be a string')
  elseif size(Param,1)~=1
    error ('Parameter must be a non-empty single row string.')
  end
  switch lower (Param)
  case 'light'
    light = lower (Value);
    if (~ischar (light))
      error ('light must be a string.')
    elseif ~(strcmp(light,'off') || strcmp(light,'on'))
      error ('light must be off | on')
    end
  case 'colormap'
    if ischar (Value)
      cmap = lower(Value);
    elseif size (Value, 2) ~= 3
      error ('colormap must be a string or have exactly three columns.')
    else
      cmap=Value;
    end
  otherwise
    error ('Unknown parameter: %s', Param)
  end
end

colormap (cmap);

% convert the number of subdivisions in number of points
subd = subd+1;

% plot the curve or surface
if (iscell (nurbs.knots))
 if (size (nurbs.knots,2) == 2) % plot a NURBS surface
  knt = nurbs.knots;
  p = nrbeval (nurbs, {linspace(knt{1}(1),knt{1}(end),subd(1)) ...
                       linspace(knt{2}(1),knt{2}(end),subd(2))});
  if (strcmp (light,'on'))
    % light surface
    surf (squeeze(p(1,:,:)), squeeze(p(2,:,:)), squeeze(p(3,:,:)),'EdgeColor','none','FaceColor',cmap);
    %shading interp;
    plotSpanLim(nurbs,subd(1));
  else 
    surf (squeeze (p(1,:,:)), squeeze (p(2,:,:)), squeeze (p(3,:,:)),'EdgeColor','none');
    %shading faceted;
    plotSpanLim(nurbs,subd(2));
  end
 elseif (size (nurbs.knots,2) == 3) % plot the boundaries of a NURBS volume
  bnd = nrbextract (nurbs);
  hold_flag = ishold;
  nrbplotx (bnd(1), subd(2:3), varargin{:});
  hold on
  %plotSpanLim(bnd(1),subd(2));
  nrbplotx (bnd(2), subd(2:3), varargin{:});
  %plotSpanLim(bnd(2),subd(2));
  nrbplotx (bnd(3), subd([1 3]), varargin{:});
  %plotSpanLim(bnd(3),subd(1));
  nrbplotx (bnd(4), subd([1 3]), varargin{:});
  %plotSpanLim(bnd(4),subd(1));
  nrbplotx (bnd(5), subd(1:2), varargin{:});
  %plotSpanLim(bnd(5),subd(1));
  nrbplotx (bnd(6), subd(1:2), varargin{:});
  %plotSpanLim(bnd(6),subd(1));
  
  if (~hold_flag)
    hold off
  end
 
 else
  error ('nrbplot: some argument is not correct')
 end
else
  % plot a NURBS curve
  p = nrbeval (nurbs, linspace (nurbs.knots(1), nurbs.knots(end), subd));

  if (any (nurbs.coefs(3,:)))
    % 3D curve
    plot3 (p(1,:), p(2,:), p(3,:)); 
    grid on;
  else
    % 2D curve
    plot (p(1,:), p(2,:));
  end
end
axis equal;

end