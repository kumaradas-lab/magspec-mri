function varargout = arrow3d(varargin)
%% Plot 3-D arrows
%
%       [hline, hhead] = arrow3d(start, stop, ang, linetype, p, n, c)
%       [hline, hhead] = arrow3d(ax, ...)
%
% This function is used to draw 3D-arrows. The size of arrow head is decided by
% the minimum distance between start and stop.
%
%
% INPUT:
%
%   ax
%       Handle to an axes graphics object. If omitted, the current axes (gca)
%       are used.
%
%   start
%       mx3 array where each row contains the starting points of the arrows.
%
%   stop
%       Array with the same size as "start" where each row contains the end
%       points of the arrows.
%
%   ang
%       Opening angle of the arrow head in degrees.
%       (Default: 30)
%
%   linetype
%       Body type of the arrow. Supported body types are 'line' or 'cylinder'.
%       (Default: 'cylinder')
%
%   p
%       p(1) is the ratio of the arrow head height to the shortest distance
%       between start and stop points, i.e., the relative length of the arrow
%       head as calculated for the shortest arrow.
%       p(2) is the ratio of the radius of the base of the arrow head to its
%       length.
%       Both values must be between 0 and 1.
%       (Default: [0.25, 0.1])
%
%   n
%       n(1) is the number of points in the base of the regular pyramid that is
%       used as the arrow head.
%       n(2) is the number of points in the base of the (elongated) regular
%       prism that is drawn as the body of the arrow if "linetype" is
%       "cylinder".
%       (Default: [20, 10])
%
%   c
%       Color value for the arrow.
%       (Default: 'k')
%
%
% OUTPUT:
%
%   hline
%       Vector with handles to the graphics object that form the bases of the
%       arrows.
%
%   hhead
%       Vector with handles to the graphics objects that form the heads of the
%       arrows.
%
%
% EXAMPLE:
%
%   t = linspace(0, 4*pi, 40);
%   x = cos(t);
%   y = sin(t);
%   z = 0.2*t;
%   p = [x.', y.', z.'];
%   p1 = p(1:end-1,:);
%   p2 = p(2:end,:);
%   hax = axes();
%   arrow3d(hax, p1, p2, 15, 'cylinder', [0.5,0.5]);
%   axis(hax, 'equal');
%   grid(hax, 'on');
%
%
% Original Function from Matlab File Exchange (license free):
%   Changshun Deng (2025). Draw 3D arrows
%   (https://www.mathworks.com/matlabcentral/fileexchange/8396-draw-3d-arrows),
%   MATLAB Central File Exchange. Retrieved June 5, 2025.
%
%   Author: Changshun Deng
%   Email: heroaq_2002@163.com
%   WEB-Log: http://waitingforme.yculblog.com
%   30/8/2005
%
%   Bug Fixed:
%   1. arrow3d([ 0 0 -1 ], [ 0 0 -2]) points the wrong way
%      Found by Pavel Grinfeld(pg@math.drexel.edu)
%      Fixed By: WaitingForMe 2006/7/24
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% check input

if nargin > 0 && ~isempty(varargin{1}) && ishghandle(varargin{1}(1), 'axes')
  hax = varargin{1}(1);
  haveAx = 1;
else
  hax = gca();
  haveAx = 0;
end

if nargin < 2 + haveAx
  error('PD:arrow3d:Nargin', ...
    'Function must be called with at least two input arguments.');
end

%% default data
start = varargin{1+haveAx};
stop = varargin{2+haveAx};

if nargin < 3 + haveAx || isempty(varargin{3+haveAx})
  ang = 30;
else
  ang = varargin{3+haveAx};
end
if nargin < 4 + haveAx || isempty(varargin{4+haveAx})
  ltype = 'cylinder';
else
  ltype = varargin{4+haveAx};
end
if nargin < 5 + haveAx || isempty(varargin{5+haveAx})
  p = [0.25, 0.1];
else
  p = varargin{5+haveAx};
end
if nargin < 6 + haveAx || isempty(varargin{6+haveAx})
  n = [20, 10];
else
  n = varargin{6+haveAx};
end
if nargin < 7 + haveAx || isempty(varargin{7+haveAx})
  c = 'k';
else
  c = varargin{7+haveAx};
end

% Check start and stop point data
if size(start, 2) ~= 3 || size(stop, 2) ~= 3
  error('PD:arrow3d:InvalidStartStop', ...
    'Start point data and stop point data must be m%c3 matrices.', char(215));
end
if any(size(start) ~= size(stop))
  error('PD:arrow3d:IncompatibleStartStop', ...
    'Start point data and stop point data must have the same size.');
end

% p(1) is the ratio of arrow head height to the distance between start and end points
% p(2) is the ratio of cylinder radial to arrow head height
% check if p is between 0 and 1
if ~all(p<=1 & p>0)
  error('PD:arrow3d:InvalidHeadSize', ...
    'Values of p must be between 0 and 1');
end

% n(1) is the patch number for arrow head
% n(2) is the patch number for cylinder
% check if n is a positive integer greater than 2
n = ceil(n);
if ~all(n>2)
  error('PD:arrow3d:InvalidNumberPatches', ...
    'Number of patches in arrow head and in cylinder must be at least 3.');
end
% direction vectors between start and end points
dvec = stop - start;
% distances between start and end points
dis = sqrt(sum(dvec.^2, 2));
% height of arrow head
hv = min(dis) * p(1);
% initial line data of 3d arrows
% init_start = zeros(size(start));
init_stop = [zeros(size(dis)), zeros(size(dis)), (dis-hv).*ones(size(dis))];
% rotate angles of the lines
cosrang = acosd(dvec(:,3)./dis);
% normal vector between arrow line and Z-axis
nvec = [-dvec(:,2), dvec(:,1), zeros(size(dis))];
% draw lines of arrows
if ~ishold(hax)
  % set hold to "on" temporarily to draw multiple arrows
  hold(hax, 'on');
  protectHold = onCleanup(@() hold(hax, 'off'));
  % view(ax,3);
end
hlines = [];
if strcmp(ltype,'line')
  for i = 1:length(dis)
    % Rotate end point
    [rx, ry, rz] = rotatedata(init_stop(i,1), init_stop(i,2), init_stop(i,3), nvec(i,:), cosrang(i), [0,0,0]);
    hlines(i) = line(hax, [start(i,1);start(i,1)+rx], [start(i,2);start(i,2)+ry], [start(i,3);start(i,3)+rz]);
  end
  hlgrd = [];
else
  for i = 1:length(dis)
    r = hv * tand(ang) .* p(2);
    [xi, yi, zi] = cylinder(hax, r.*[1,1], n(2));
    zi = zi .* (dis(i)-hv);
    % escape the error if the arrow is in z-direction
    % if the arrow is in z-direction then the nvector result zeros to
    % make a error!
    % Fix this bug! 2006/07/24
    % Thanks to Pavel Grinfeld(pg@math.drexel.edu)
    if all(nvec(i,:) == 0)
      nvec(i,:) = [0, 1, 0];
    end
    [rx, ry, rz] = rotatedata(xi ,yi, zi, nvec(i,:), cosrang(i), [0,0,0]);
    cx = start(i,1) + rx;
    cy = start(i,2) + ry;
    cz = start(i,3) + rz;
    hlines(i) = surf(hax, cx, cy, cz, 'EdgeColor', 'none', 'FaceColor', c);
    hlgrd(i) = patch(hax, cx(1,:), cy(1,:), cz(1,:), c, 'EdgeColor', 'none');
  end
end

% initial arrow head data of 3d arrows
hheads = [];
hhgrd = [];
pv = dis-hv;
% draw heads of arrows
for i = 1:length(dis)
  % initial taper data
  [xi, yi, zi] = cylinder(hax, [tand(ang), 0], n(1));
  xi = xi*hv;
  yi = yi*hv;
  zi = zi*hv+pv(i);
  % rotate taper
  [rx, ry, rz] = rotatedata(xi, yi, zi, nvec(i,:), cosrang(i), [0,0,0]);
  cx = start(i,1) + rx;
  cy = start(i,2) + ry;
  cz = start(i,3) + rz;
  hheads(i) = surf(hax, cx, cy, cz, 'EdgeColor', 'none', 'FaceColor', c);
  % draw underside of taper
  hhgrd(i) = patch(hax, cx(1,:), cy(1,:), cz(1,:), c, 'EdgeColor', 'k');
end

if nargout>1
  varargout{2} = [hheads; hhgrd];
end
if nargout>0
  varargout{1} = [hlines; hlgrd];
end

end


function [newx, newy, newz] = rotatedata(xdata, ydata, zdata, azel, alpha, origin)
%% ROTATEDATA rotate data about specified origin and direction.
%
%   ROTATEDATA(Xdata,Ydata,Zdata,[THETA PHI],ALPHA,ORIGIN) rotates the objects with handles H
%   through angle ALPHA about an axis described by the 2-element
%   direction vector [THETA PHI] (spherical coordinates).
%   All the angles are in degrees.  The handles in H must be children
%   of the same axes.
%
%   THETA is the angle in the xy plane counterclockwise from the
%   positive x axis.  PHI is the elevation of the direction vector
%   from the xy plane (see also SPH2CART).  Positive ALPHA is defined
%   as the righthand-rule angle about the direction vector as it
%   extends from the origin.
%
%   ROTATEDATA(Xdata,Ydata,Zdata,[X Y Z],ALPHA,ORIGIN) rotates the objects about the direction
%   vector [X Y Z] (Cartesian coordinates). The direction vector
%   is the vector from the center of the plot box to (X,Y,Z).
%
%   See also SPH2CART, CART2SPH.
%
%   Modified by ChangShun Deng
%   Email: heroaq_2002@163.com
%   Web-Log: http://waitingforme.yculblog.com
%   2005/3/4
%
%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 5.17 $  $Date: 2002/06/05 20:05:16 $
%%%%

if nargin < 6
  error('Not enough input arguments! Type ''help rotatedata'' to get some help!')
end

if numel(azel) == 2  % theta, phi
    theta = azel(1);
    phi = azel(2);
    u = [cosd(phi)*cosd(theta); cosd(phi)*sind(theta); sind(phi)];
elseif numel(azel) == 3  % direction vector
    u = azel(:)/norm(azel);
end

cosa = cosd(alpha);
sina = sind(alpha);
vera = 1 - cosa;
x = u(1);
y = u(2);
z = u(3);

rot = [cosa+x^2*vera,   x*y*vera-z*sina, x*z*vera+y*sina; ...
       x*y*vera+z*sina, cosa+y^2*vera,   y*z*vera-x*sina; ...
       x*z*vera-y*sina, y*z*vera+x*sina, cosa+z^2*vera]';
[m, n] = size(xdata);
% if isempty(z)
%   z = zeros(size(x));
% end
newxyz = [xdata(:)-origin(1), ydata(:)-origin(2), zdata(:)-origin(3)];
newxyz = newxyz*rot;
newx = origin(1) + reshape(newxyz(:,1),m,n);
newy = origin(2) + reshape(newxyz(:,2),m,n);
newz = origin(3) + reshape(newxyz(:,3),m,n);

end
