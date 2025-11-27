function varargout = plotMinMaxInterp(varargin)
%% Faster plotting of large timeseries
%
%   plotMinMaxInterp(hAxes, x, y, PropertyName, PropertyValue)
%   plotMinMaxInterp(y)
%   plotMinMaxInterp(x, y)
%   plotMinMaxInterp(x, y, PropertyName, PropertyValue)
%   h = plotMinMaxInterp(...)
%
% This function is a minimal wrapper for the Matlab plot function with automatic
% downsampling at low zoom factors and cropping at high zoom factors for faster
% zoom and pan.
% It works in subplots and can be combined with normals plots of any kind in the
% same axes.
%
% The amount of downsampling can be changed by adjusting the value of MAXPOINTS.
% Lines with less than MAXPOINTS are interpolated to approx. that number of
% points.
% Plotting of multiple lines at once is allowed up to MAXHANDLES. Adjust as
% necessary.
%
% Note that plotMinMaxInterp uses the 'Tag' and 'UserData' properties of the
% plotted lines.
%
% Derived from the function "jplot" from the Matlab File Exchange:
% https://www.mathworks.com/matlabcentral/fileexchange/42191-jplot
% v1.0 Jake Reimer 6/11/2013
% v1.1 Jake Reimer 6/12/2013 Takes initial xL as [min(x) max(x)] rather than
%                            initial xlims so that downsampling works even if
%                            you're zoomed into the axis when you jplot.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------


%% constants (change as needed)
MAXPOINTS = 2^14;
MAXHANDLES = 1000;


%% check input arguments
hAxes = [];
if nargin == 0
  return;
elseif nargin == 1
  if isa(varargin{1}, 'handle')
    return;
  else
    y = varargin{1};
    x = 1:length(y);
    mods = {};
  end
elseif nargin == 2
  if isa(varargin{1}, 'handle')
    if ~isvalid(varargin{1})
      return;
    end
    hAxes = varargin{1};
    y = varargin{2};
    x = 1:length(y);
    mods = {};
  else
    x = varargin{1};
    y = varargin{2};
    mods = {};
  end
elseif nargin == 3
  if isa(varargin{1}, 'handle')
    if ~isvalid(varargin{1})
      return;
    end
    hAxes = varargin{1};
    x = varargin{2};
    y = varargin{3};
    mods = {};
  else
    x = varargin{1};
    y = varargin{2};
    mods = varargin(3:end);
  end
elseif nargin > 3
  if isa(varargin{1}, 'handle')
    if ~isvalid(varargin{1})
      return;
    end
    hAxes = varargin{1};
    x = varargin{2};
    y = varargin{3};
    mods = varargin(4:end);
  else
    x = varargin{1};
    y = varargin{2};
    mods = varargin(3:end);
  end
end
if isempty(hAxes)
  hAxes = get(gcf, 'CurrentAxes');
end

if isvector(y)
  numLines = 1;
else
  numLines = size(y, 2);
end
if numLines > MAXHANDLES
  error('PD:plotMinMaxInterp:TooManyLines', ...
    'Number of lines to plot (%d) exceeds maximum (%d).', numLines, MAXHANDLES);
end

if isvector(y)
  y = reshape(y, [], 1);
end
if isvector(x)
  x = repmat(reshape(x, [], 1), 1, size(y, 2));
end


%% set zoom and pan callbacks
hZoom = zoom(gcf);
hPan = pan(gcf);
set(hZoom, 'ActionPostCallback', @zoomCallback);
set(hPan , 'ActionPostCallback', @zoomCallback);


%% prepare data for plot
yplot = NaN(MAXPOINTS, size(y, 2));
xplot = NaN(MAXPOINTS, size(y, 2));
raw(1:size(y, 2)) = struct('len', [], 'x', [], 'y', [], 'xL', [], 'maxPoints', []);
for n = 1:size(y, 2)
  raw(n).len = size(x, 1);
  raw(n).x = x(:, n);
  raw(n).y = y(:, n);
  raw(n).xL = [min(raw(n).x), max(raw(n).x)];
  raw(n).maxPoints = MAXPOINTS;

  [xplot(:,n), yplot(:,n)] = xyPlot(raw(n).x, raw(n).y, MAXPOINTS);
end


%% plot data
if ~ishold(hAxes)
  cla(hAxes);
end
hLines(1:size(y, 2)) = matlab.graphics.GraphicsPlaceholder;
for n = 1:size(y, 2)
  hLines(n) = line(xplot(:,n), yplot(:,n), 'Parent', hAxes, mods{:});
  set(hLines(n), 'UserData', raw(n), 'Tag', 'dec');
end


%% (optional) output argument
if nargout > 0
  varargout{1} = hLines;
end

end


function [xplot, yplot] = xyPlot(xSel, ySel, maxPoints)
%% Interpolate the function to get between 4 and maxPoints/4 nodes
% (at the current zoom level)

% Interp 2^2 minimum to find proper min max values
os2 = max(min(floor(log2(maxPoints/numel(xSel))), 12), 2);

xplot = reshape(interpn(xSel, os2, 'linear'), [], 1);
yplot = reshape(interp1(xSel, ySel, xplot, 'spline'), [], 1);

if numel(xSel) > maxPoints/4
  % plot a block of vertical line segments
  % (to indicate that some information isn't displayed at the current zoom level)
  len = numel(xplot);
  dec = ceil(len/(maxPoints/2));  % find segment length to get half maxPoints segments
  xstart = reshape(xplot(1:dec:floor(len/dec)*dec), 1, []);
  xplot = reshape([xstart; xstart] + max(diff(xstart(1:2)), 0)/2, floor(len/dec)*2, 1);
  temp = reshape(yplot(1:floor(len/dec)*dec), dec, floor(len/dec));
  % find min and max per segment
  yplot = reshape([max(temp, [], 1); min(temp, [], 1)], floor(len/dec)*2, 1);
end

xplot(end:maxPoints) = NaN;
yplot(end:maxPoints) = NaN;

end


function zoomCallback(~, ev)
%% Update plot on zoom or pan

if strcmp(get(gcf, 'selectiontype'), 'normal')
  xL = get(ev.Axes, 'xlim');
else
  xL = [-Inf, Inf];
end

h = findobj(ev.Axes, 'Tag', 'dec');

for n = 1:length(h)
  raw = get(h(n), 'UserData');
  ind = [max(find(raw.x>=xL(1)-diff(xL)/2, 1, 'first'), 1), ...
    min(find(raw.x<=xL(2)+diff(xL)/2, 1, 'last'), numel(raw.x))];  % expand half width
  ind = [max(ind(1)-3, 1), min(ind(2)+3, numel(raw.x))];  % expand -/+ 2 samples
  xSel = raw.x(ind(1):ind(2));
  ySel = raw.y(ind(1):ind(2));
  maxPoints = raw.maxPoints;
  [xplot, yplot] = xyPlot(xSel, ySel, maxPoints);
  set(h(n), 'XData', xplot, 'YData', yplot);
end

drawnow();

end
