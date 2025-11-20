function hIso = DrawLocalIsoSurface(this, volume, data, datanormals, value, oldiso)
%% Draw or refresh iso surface

pushset(this.hFigure, 'Pointer', 'watch');

fv = isosurface(volume{:}, data, value);

clim = get(this.hAxes, 'CLim');
cmap = get(this.hFigure, 'ColorMap');
clen = clim(2)-clim(1);
idx = fix((value-clim(1))*length(cmap)/clen)+1;

if nargin > 5
  try
    set(oldiso, fv, 'FaceColor', cmap(max(1,min(end,idx)),:));
  catch
    set(oldiso, fv, 'FaceColor', 'none');
  end
  hIso = oldiso;
  cap = getappdata(hIso, 'isosurfacecap');
  if ~isempty(cap)
    this.DrawLocalIsoCaps(hIso, cap);
  end
else
  hIso = patch(fv, 'EdgeColor', 'none', 'FaceColor', cmap(max(1,min(end,idx)),:), ...
    'Tag', 'sliceomaticisosurface', 'Parent', this.hAxes);
  switch get(this, 'deflight')
    case 'flat'
      set(hIso, 'FaceLighting', 'flat');
    case 'smooth'
      set(hIso, 'FaceLighting', 'gouraud');
  end
  setappdata(hIso, 'isosurfacecap', []);
end

setappdata(hIso, 'isosurfacevalue', value);
setappdata(hIso, 'isosurfacedata', data);

% disable lighting for this patch before reducing vertices
fl = get(hIso, 'FaceLighting');
el = get(hIso, 'EdgeLighting');
set(hIso, 'FaceLighting', 'none', 'EdgeLighting', 'none', 'VertexNormals', []);

% reduce vertices (does not update VertexNormals)
reducepatch(hIso, 10000);

% calculate normals with reduced number of vertices and enable lighting again
isonormals(volume{:}, datanormals, hIso);
set(hIso, 'FaceLighting', fl, 'EdgeLighting', el);

popset(this.hFigure, 'Pointer');

end
