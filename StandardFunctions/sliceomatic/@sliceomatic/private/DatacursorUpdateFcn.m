function dcText = DatacursorUpdateFcn(~, eventObj)
%% Function for text at Datacursor

%% default output
dcText = '';

%% check object
dcTarget = get(eventObj, 'Target');
hAxes = ancestor(dcTarget, 'axes');
if ~strcmpi(get(hAxes, 'Tag'), 'MainAxes')
  return
end

hParent = get(hAxes, 'Parent');
d = getappdata(hParent, 'sliceomatic');
this = d.axmain_obj;
dcPosition = get(eventObj, 'Position');


%% get value from position
if ishghandle(dcTarget, 'surface') % slice
  % FIXME: probably better to interpolate data
  if isnan(this.xmesh)
    xi = round(dcPosition(1));
  else
    dist = abs(this.xmesh - dcPosition(1));
    xi = find(min(dist) == dist, 1);
  end
  if isnan(this.ymesh)
    yi = round(dcPosition(2));
  else
    dist = abs(this.ymesh - dcPosition(2));
    yi = find(min(dist) == dist, 1);
  end
  if isnan(this.zmesh)
    zi = round(dcPosition(3));
  else
    dist = abs(this.zmesh - dcPosition(3));
    zi = find(min(dist) == dist, 1);
  end
  
  % get value at position
  value = d.data(yi,xi,zi);
elseif isappdata(dcTarget, 'isosurfacevalue') % iso-surface
  value = getappdata(dcTarget, 'isosurfacevalue');
else
  return
end
dcText = sprintf('X: %.3g\nY: %.3g\nZ: %.3g\nV: %.3g', dcPosition, value);

end
