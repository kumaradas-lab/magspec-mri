function ArrowDragFinish(this)
%% Finish dragging an arrow. Finalize graphics objects

arrows = findall(this.hParent, 'Tag', 'sliceomaticarrow');

popset(this.draggedArrow, 'FaceColor');
popset(this.draggedArrow, 'FaceAlpha');

popset(arrows, 'FaceColor');
popset(arrows, 'FaceAlpha');

ss = getappdata(this.draggedArrow, 'arrowslice');
if isempty(ss)
  ss = getappdata(this.draggedArrow, 'arrowiso');
end

% These pushes are junk which will be undone when all slices or
% isosurfs are reset below.
pushset(ss, 'FaceAlpha', 1);
pushset(ss, 'EdgeColor', 'k');

slices = this.GetAllSlices();

if ~isempty(slices)
  popset(slices, 'FaceAlpha');
  popset(slices, 'EdgeColor');
end

isosurfs = this.GetAllIsos();

if ~isempty(isosurfs)
  popset(isosurfs, 'FaceAlpha');
  popset(isosurfs, 'EdgeColor');
end

data = get(this, 'data');
smooth = get(this, 'smooth');

if isnan(this.xmesh)
  for iIso = 1:length(isosurfs)
    cap = getappdata(isosurfs(iIso), 'isosurfacecap');
    if ~isempty(cap)
      popset(cap, 'Visible');
      this.DrawLocalIsoCaps(isosurfs(iIso), cap);
    end
    if getappdata(isosurfs(iIso), 'reduced')
      setappdata(isosurfs(iIso), 'reduced', 0);
      this.DrawLocalIsoSurface({}, data, smooth, ...
        getappdata(isosurfs(iIso), 'isosurfacevalue'), ...
        isosurfs(iIso));
    end
  end
else
  for iIso = 1:length(isosurfs)
    cap = getappdata(isosurfs(iIso), 'isosurfacecap');
    if ~isempty(cap)
      popset(cap, 'Visible');
      this.DrawLocalIsoCaps(isosurfs(iIso), cap);
    end
    if getappdata(isosurfs(iIso), 'reduced')
      setappdata(isosurfs(iIso), 'reduced', 0);
      realvolume = {this.xmesh, this.ymesh, this.zmesh};
      this.DrawLocalIsoSurface(realvolume, data, smooth, ...
        getappdata(isosurfs(iIso), 'isosurfacevalue'), ...
        isosurfs(iIso));
    end
  end
end

popset(this.hFigure, 'WindowButtonUpFcn');
popset(this.hFigure, 'WindowButtonMotionFcn');

showarrowtip(this.hParent, []);

% Make sure whatever buttonupfcn on the figure is run now to "turn
% off" whatever was going on before we got our callback on the
% arrow.

buf = get(this.hFigure, 'WindowButtonUpFcn');
checkEvalBuffer(buf)

% re-set figure mode
tempfiguremode(this.hFigure, 'reset');

set(this.hFigure, 'Pointer', 'arrow');

end
