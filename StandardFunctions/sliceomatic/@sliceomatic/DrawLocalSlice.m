function hSurface = DrawLocalSlice(this, data, X, Y, Z, oldslice)
%% Slice Management
% Uses specialized slicomatic slices, not slices created with the SLICE command.
%
% -----------------------------------------------------------------------
% (C) Copyright 2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% -----------------------------------------------------------------------

if ~isempty(X)
  [xdata, ydata, zdata, cdata] = getDataSlice(this, data, X, 'x');
  slice_type = 'X';
elseif ~isempty(Y)
  [xdata, ydata, zdata, cdata] = getDataSlice(this, data, Y, 'y');
  slice_type = 'Y';
elseif ~isempty(Z)
  [xdata, ydata, zdata, cdata] = getDataSlice(this, data, Z, 'z');
  slice_type = 'Z';
else
  % FIXME: Can we ever reach here
  error('Incorrect syntax for sliceomatic.DrawLocalSlice.');
end

cdata = squeeze(cdata);
xdata = squeeze(xdata);
ydata = squeeze(ydata);
zdata = squeeze(zdata);

if nargin > 5
  % Recycle the old slice
  set(oldslice, 'CData', cdata, 'AlphaData', cdata, 'XData', xdata, ...
    'YData', ydata, 'ZData', zdata);
  hSurface = oldslice;
  %delete(news);
  if propcheck(hSurface, 'FaceColor', 'texturemap') && ...
      ((all(diff(xdata(:))==0) && strcmpi(get(this.hAxes, 'YScale'), 'linear') && ...
      strcmpi(get(this.hAxes, 'ZScale'), 'linear')) ...
      || (all(diff(ydata(:))==0) && strcmpi(get(this.hAxes, 'YScale'), 'linear') && ...
      strcmpi(get(this.hAxes, 'XScale'), 'linear')) ...
      || (all(diff(zdata(:))==0) && strcmpi(get(this.hAxes, 'XScale'), 'linear') && ...
      strcmpi(get(this.hAxes, 'YScale'), 'linear')))
    textureizeslice(hSurface, 'on');
  end
  setappdata(hSurface, 'slicetype', slice_type);
else
  % setup the alphadata
  new_surface = surface('CData', cdata, 'XData', xdata, ...
    'YData', ydata, 'ZData', zdata, 'Parent', this.hAxes);
  set(new_surface, 'AlphaData', cdata, 'AlphaDataMapping', 'scaled', 'Tag', 'sliceomaticslice', ...
    'FaceLighting', 'none', ...
    'UIContextMenu', this.sliceContextMenu);
  hSurface = new_surface;
  setappdata(hSurface, 'slicetype', slice_type);
  switch get(this, 'defcolor')
    case 'faceted'
      set(hSurface, 'FaceColor', 'flat', 'EdgeColor', 'k');
    case 'flat'
      set(hSurface, 'FaceColor', 'flat', 'EdgeColor', 'none');
    case 'interp'
      set(hSurface, 'FaceColor', 'interp', 'EdgeColor', 'none');
    case 'texture'
      set(hSurface, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
  end
  switch get(this, 'defalpha')
    case 'none'
      set(hSurface, 'FaceAlpha', 1);
    case 'flat'
      set(hSurface, 'FaceAlpha', 'flat');
    case 'interp'
      set(hSurface, 'FaceAlpha', 'interp');
    case 'texture'
      set(hSurface, 'FaceAlpha', 'texturemap');
  end
end

contour = getappdata(hSurface, 'contour');
if ~isempty(contour)
  try
    levels = getappdata(hSurface, 'contourlevels');
    if ~isempty(levels)
      this.DrawLocalContour(hSurface, contour, levels);
    else
      this.DrawLocalContour(hSurface, contour);
    end
  catch
    this.DrawLocalContour(hSurface, contour);
  end
end

end


function [xdata, ydata, zdata, cdata] = getDataSlice(this, data, coor, type)
%% Get data for selected slice

meshField = [type 'mesh'];
ds = size(data);
iType = type-119; % 'x' -> 1, 'y' -> 2, 'z' -> 3

if isnan(this.(meshField))
  gridVecs = {1:ds(2), 1:ds(1), 1:ds(3)};
  idx = round(coor);
  gridVecs{iType} = idx;
  [xdata, ydata, zdata] = meshgrid(gridVecs{:});
  if idx > 0 && idx <= gridVecs{iType}(end)
    switch type
      case 'x'
        cdata = reshape(data(:,idx,:), ds(1), ds(3));
      case 'y'
        cdata = reshape(data(idx,:,:), ds(2), ds(3));
      case 'z'
        cdata = reshape(data(:,:,idx), ds(1), ds(2));
    end
  else
    cdata = nan(size(xdata));
  end
else
  diffCoor = abs(this.(meshField) - coor);
  slice_number = find(diffCoor == min(diffCoor) & ...
    diffCoor <= abs(this.(meshField)(1) - this.(meshField)(2))/2, 1);
  if isequal(get(this, 'zdir'), 'reverse')
    slice_number = length(this.(meshField)) - slice_number + 1;
  end
  switch type
    case 'x'
      [xdata, ydata, zdata] = meshgrid(this.(meshField)(slice_number), this.ymesh, this.zmesh);
    case 'y'
      [xdata, ydata, zdata] = meshgrid(this.xmesh, this.(meshField)(slice_number), this.zmesh);
    case 'z'
      [xdata, ydata, zdata] = meshgrid(this.xmesh, this.ymesh, this.(meshField)(slice_number));
  end
  iTypeSize = mod(2-(1:3),3)+1;
  if ~isempty(slice_number) && slice_number(1) > 0 && slice_number(1) <= ds(iTypeSize(iType))
    switch type
      case 'x'
        cdata = reshape(data(:,slice_number(1),:), ds(1), ds(3));
      case 'y'
        cdata = reshape(data(slice_number(1),:,:), ds(2), ds(3));
      case 'z'
        cdata = reshape(data(:,:,slice_number(1)), ds(1), ds(2));
    end
  else
    cdata = nan(size(xdata));
  end
end

end
