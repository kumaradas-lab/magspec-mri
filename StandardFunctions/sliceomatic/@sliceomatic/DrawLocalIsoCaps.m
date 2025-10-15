function hIsoCap = DrawLocalIsoCaps(this, hIsoSurface, oldIsoCap)
%% Draw or refresh caps for iso surfaces

% Get relevant info from the isosurface.
if nargin < 3 || ~strcmp(get(oldIsoCap, 'Visible'), 'off')
  data = getappdata(hIsoSurface, 'isosurfacedata');
  if isnan(this.xmesh)
    caps = isocaps(data, getappdata(hIsoSurface, 'isosurfacevalue'));
  else
    caps = isocaps(this.xmesh, this.ymesh, this.zmesh, data, getappdata(hIsoSurface, 'isosurfacevalue'));
  end
end

if nargin > 2
  if ~strcmp(get(oldIsoCap, 'Visible'), 'off')
    % FIXME: When do we refrresh the data if visible is currently off?
    set(oldIsoCap, caps);
  end
  hIsoCap = oldIsoCap;
else
  hIsoCap = patch(caps, 'EdgeColor', 'none', 'FaceColor', 'flat', ...
    'FaceLighting', 'none', ...
    'Tag', 'sliceomaticisocap', 'Parent', this.hAxes);

  setappdata(hIsoCap, 'isosurface', hIsoSurface);
  setappdata(hIsoSurface, 'isosurfacecap', hIsoCap);

  switch get(this, 'defcolor')
    case 'faceted'
      set(hIsoCap, 'FaceColor', 'flat', 'EdgeColor', 'black');
    case 'flat'
      set(hIsoCap, 'FaceColor', 'flat', 'EdgeColor', 'none');
    case 'interp'
      set(hIsoCap, 'FaceColor', 'interp', 'EdgeColor', 'none');
    case 'texture'
      set(hIsoCap, 'FaceColor', 'flat', 'EdgeColor', 'none');
    case 'none'
      set(hIsoCap, 'FaceColor', 'none', 'EdgeColor', 'none');
  end
  switch get(this, 'defalpha')
    case 'none'
      set(hIsoCap, 'FaceAlpha', 1);
    case 'flat'
      set(hIsoCap, 'FaceAlpha', 'flat');
    case 'interp'
      set(hIsoCap, 'FaceAlpha', 'interp');
    case 'texture'
      set(hIsoCap, 'FaceAlpha', 'flat');
  end
end

end
