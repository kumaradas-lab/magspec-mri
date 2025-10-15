function hLines = sheerPlot(varargin)
%% Plot staggered lines on sheered axes
%
%   hLines = sheerPlot(Y)
%   hLines = sheerPlot(X, Y)
%   hLines = sheerPlot(X, Y, Z)
%   hLines = sheerPlot(hAxes, ...)
%   hLines = sheerPlot(..., 'prop', val, ...)
%
%
% INPUT:
%
%   X
%         Column vector with the X data (i.e. the axis that points to the right)
%         Default: 1:size(Y, 1)
%
%   Y
%         Matrix with the Y values of the sheered lines. Each row of the matrix
%         contains the data of one line. The number of columns must match X.
%
%   Z
%         Vector with values along the sheered axis. It must contain one value
%         for each line. Default: 1:size(Y, 2)
%
%   hAxes
%         Handle to the axes that is used to plot the sheered lines. If omitted,
%         a new axes object is created.
%
%   'prop', val
%         Property value pairs. These properties are passed on when creating the
%         line objects apart from the following that have special meaning:
%     'OffsetX'
%           Scalar with offset in x direction (to the right) between the lines
%           in parts of the total line length. Default: 0.1
%     'OffsetY'
%           Scalar with offset in y direction (to the top) between the lines in
%           parts of the maximum range of the lines. Default: 1
%     'XLabel'
%           Label on the x-axis (horizontal). Default: ''
%     'ZLabel'
%           Label on the sheered axis. Default: ''
%
% OUTPUT:
%
%   hLines
%         Handles to the created lines.
%
%
% PROGRAMMING NOTES:
%
%   The axes are drawn manually. So zooming, panning and rotating the underlying
%   axes might lead to unexpected results. Also use functions carefully that
%   manipulate the underlying axes. Instead, set all necessary properties when
%   calling "sheerPlot".
%
% ------------------------------------------------------------------------------
% (C) Copyright 2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% input check
% axes handle
if nargin > 0
  if ishghandle(varargin{1}, 'axes')
    hAxes = varargin{1};
    cla(hAxes);
    varargin(1) = [];
  else
    hAxes = gca;
  end
end

% input X, Y, and Z
if isempty(varargin) || ~isnumeric(varargin{1})
  haveInputData = false;
else
  haveInputData = true;
  if numel(varargin) > 1 && isnumeric(varargin{2})
    X = varargin{1};
    Y = varargin{2};
    if numel(varargin) > 2 && isnumeric(varargin{3})
      Z = varargin{3};
      varargin(1:3) = [];
    else
      Z = 1:size(Y,2);
      varargin(1:2) = [];
    end
  else
    Y = varargin{1};
    X = 1:size(Y, 1);
    Z = 1:size(Y, 2);
    varargin(1) = [];
  end
end

if ~haveInputData
  % check if data is stored in axes
  if isappdata(hAxes, 'XData') && isappdata(hAxes, 'YData') && isappdata(hAxes, 'ZData')
    X = getappdata(hAxes, 'XData');
    Y = getappdata(hAxes, 'YData');
    Z = getappdata(hAxes, 'ZData');
  else
    error('PD:sheerPlot:NoData', 'sheerPlot needs at least one numeric input.');
  end
end

if size(X, 1) ~= size(Y, 1)
  error('PD:sheerPlot:SizeMismatchRows', ...
    'Input X and Y must have the same number of rows');
end

if size(X, 2) == 1
  X = repmat(X, [1, size(Y, 2)]);
end
if size(X, 2) ~= size(Y, 2)
  error('PD:sheerPlot:SizeMismatchColumnsX', ...
    'Input X must be a column vector or have the same number of columns as input Y.');
end

if numel(Z) ~= size(Y, 2)
  error('PD:sheerPlot:SizeMismatchColumnsZ', ...
    'Input Z must have one element for each column in input Y.');
end


%% special sheerPlot property value pairs

% default values
offsetX = 0.1;
offsetY = 1;
xLabel = '';
zLabel = '';

% restore values from previous axes if possible
if isappdata(hAxes, 'OffsetX')
  offsetX = getappdata(hAxes, 'OffsetX');
end
if isappdata(hAxes, 'OffsetY')
  offsetY = getappdata(hAxes, 'OffsetY');
end
if isappdata(hAxes, 'XLabel')
  xLabel = getappdata(hAxes, 'XLabel');
end
if isappdata(hAxes, 'ZLabel')
  zLabel = getappdata(hAxes, 'ZLabel');
end

% get special properties from input parameters
if numel(varargin) > 0
  if mod(numel(varargin), 2)
    error('PD:sheerPlot:PropValMisMatch', ...
      'Input must end with property value pairs.');
  end
  iProp = 0;
  while true
    if numel(varargin) < 2*iProp + 2
      break;
    end
    if ~ischar(varargin{2*iProp+1})
      error('PD:sheerPlot:InvalidProp', 'Properties must be strings.');
    end
    switch lower(varargin{2*iProp+1})
      case 'offsetx'
        if ~isnumeric(varargin{2*iProp+2}) || ~isscalar(varargin{2*iProp+2})
          error('PD:sheerPlot:InvalidOffsetX', '"OffsetX" must be a numeric scalar.');
        end
        offsetX = varargin{2*iProp+2};
        varargin(2*iProp+(1:2)) = [];
        continue;
      case 'offsety'
        if ~isnumeric(varargin{2*iProp+2}) || ~isscalar(varargin{2*iProp+2})
          error('PD:sheerPlot:InvalidOffsetY', '"OffsetY" must be a numeric scalar.');
        end
        offsetY = varargin{2*iProp+2};
        varargin(2*iProp+(1:2)) = [];
        continue;
      case 'xlabel'
        if ~ischar(varargin{2*iProp+2})
          error('PD:sheerPlot:InvalidXLabel', '"XLabel" must be a string.');
        end
        xLabel = varargin{2*iProp+2};
        varargin(2*iProp+(1:2)) = [];
        continue;
      case 'zlabel'
        if ~ischar(varargin{2*iProp+2})
          error('PD:sheerPlot:InvalidZLabel', '"ZLabel" must be a string.');
        end
        zLabel = varargin{2*iProp+2};
        varargin(2*iProp+(1:2)) = [];
        continue;
    end
    iProp = iProp + 1;
  end
end


%% store data and settings in axes
setappdata(hAxes, 'XData', X);
setappdata(hAxes, 'YData', Y);
setappdata(hAxes, 'ZData', Z);

setappdata(hAxes, 'OffsetX', offsetX);
setappdata(hAxes, 'OffsetY', offsetY);
setappdata(hAxes, 'XLabel', xLabel);
setappdata(hAxes, 'ZLabel', zLabel);


%% plot lines
Y_orig = Y;
offsetX = (max(X(:)) - min(X(:))) * offsetX;
offsetY = (max(Y(:)) - min(Y(:))) * offsetY;
maxX = max(X(:));
X = bsxfun(@plus, X, cumsum(offsetX*ones(1, size(X,2)), 2) - offsetX);
Y = bsxfun(@plus, Y, cumsum(offsetY*ones(1, size(Y,2)), 2) - offsetY);

hLines = line(hAxes, X, Y, varargin{:});

axesColor = get(groot, 'DefaultAxesColor');
xColor = get(groot, 'DefaultAxesXColor');
set(hAxes, 'Box', 'off', 'YTick', [], 'YColor', 'none', 'XColor', 'none', 'Color', 'none');


%% manually plot a grid
gridColor = get(hAxes, 'GridColor');
gridAlpha = get(hAxes, 'GridAlpha');

% y grid
xTicks = get(hAxes, 'XTick');
yLims = get(hAxes, 'YLim');
dz = diff(Z(:));
if ~isempty(dz) && all(dz-dz(1)<100*eps(dz(1)))
  yLims(1) = 0*offsetY + floor(min(Y_orig(:))/offsetY-0.01)*offsetY;
  yLims(2) = ceil((offsetY*(size(Y,2)-1)+max(Y_orig(:)))/offsetY)*offsetY;
  set(hAxes, 'YLim', yLims);
end
iMaxX = find(xTicks < maxX, 1, 'last');
if numel(xTicks) > iMaxX, xTicks(iMaxX+2:end) = []; end
set(hAxes, 'XTick', xTicks);
set(hAxes, 'XTickLabelMode', 'auto');
xValues = [xTicks + yLims(1)*offsetX/offsetY; xTicks + yLims(2)*offsetX/offsetY];
yValues = repmat(yLims.', [1, size(xValues, 2)]);

hYGrid = surface(hAxes, xValues, yValues, zeros(size(xValues)), ...
  'FaceColor', 'none', 'EdgeColor', gridColor, 'EdgeAlpha', gridAlpha);
lowerX = xValues(1,:);
upperXLim = xValues(2,[1,end]);

% x grid
xValues = bsxfun(@plus, reshape(xTicks([1,end]), [], 1), cumsum(offsetX*ones(1, size(X,2)), 2) - offsetX);
yValues = cumsum(offsetY*ones(1, size(Y,2)), 2) - offsetY;
hXGrid = surface(hAxes, xValues, yValues([1, 1], :), zeros(size(xValues)), ...
  'FaceColor', 'none', 'EdgeColor', gridColor, 'EdgeAlpha', gridAlpha);

% "axes" background
axesVertices = [lowerX(1), yLims(1); lowerX(end), yLims(1); ...
  upperXLim(end), yLims(2); upperXLim(1), yLims(2); NaN, NaN];
hBackground = patch(hAxes, 'Vertices', axesVertices, 'Faces', [1 2 3 4], ...
  'FaceColor', axesColor, 'EdgeColor', 'none');
hYAxis = patch(hAxes, 'Vertices', axesVertices, 'Faces', [4 1 2 3], ...
  'FaceColor', 'none', 'EdgeColor', xColor);

set(hAxes, 'XLim', [min([upperXLim(1), lowerX(1), xTicks(1)]), max(upperXLim(2), lowerX(end))]);


%% order children such that background and grid are plotted behind the lines
hKids = get(hAxes, 'Children');
hKids(hKids == hBackground) = [];
hKids(hKids == hXGrid) = [];
hKids(hKids == hYGrid) = [];
hKids(hKids == hYAxis) = [];

set(hAxes, 'Children', [hKids(end:-1:1); hYAxis; hXGrid; hYGrid; hBackground]);


%% manually add axes labels

% x ticks
hXLabel = get(hAxes, 'XLabel');
xLabelPos = get(hXLabel, 'Position');
yLims = get(hAxes, 'YLim');
text(hAxes, lowerX, repmat(yLims(1)-(yLims(1)-xLabelPos(2))/2, [1, size(xTicks,2)]), ...
  get(hAxes, 'XTickLabel'), ...
  'HorizontalAlignment', 'center', ...
  'VerticalAlignment', 'middle');

if ~isempty(xLabel)
  text(hAxes, mean(lowerX), xLabelPos(2), xLabel, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', get(hXLabel, 'FontSize'));
end

% z ticks
set(hAxes, 'YTick', 0, 'YTickLabel', sprintf('%g', Z(1)));
hYLabel = get(hAxes, 'YLabel');
yLabelPos = get(hYLabel, 'Position');
set(hAxes, 'YTick', []);
xLims = get(hAxes, 'XLim');

zTicksX = xValues(2,:);
zTicksY = yValues;
text(hAxes, zTicksX + 0.5*(xLims(1) - yLabelPos(1)), zTicksY, ...
  arrayfun(@(x) sprintf('%g', x), Z, 'UniformOutput', false));

set(hAxes, 'OuterPosition', [0, 0.1, 1, 0.8]);

if ~isempty(zLabel)
  % position
  xLims = get(hAxes, 'XLim');
  % angle
  dx = axesVertices(3,1) - axesVertices(2,1);
  dy = axesVertices(3,2) - axesVertices(2,2);
  dar = get(hAxes, 'DataAspectRatio');
  pbar = get(hAxes, 'PlotBoxAspectRatio');
  text(hAxes, mean(zTicksX) + 1.5*(xLims(1) - yLabelPos(1)), mean(zTicksY), zLabel, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'Rotation', atan2d(dy/(dar(2)/pbar(2)), dx/(dar(1)/pbar(1))), ...
    'FontSize', get(hYLabel, 'FontSize'));
end


%% set resize function
% FIXME: Can we adapt the layout instead of calling the function on each resize?
hParent = get(hAxes, 'Parent');
set(hParent, 'ResizeFcn', @(h,e) sheerPlot(hAxes, varargin{:}));


end
