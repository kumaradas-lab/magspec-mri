function ArrowDragPrepare(this)
%% Prepare graphics objects for dragging an arrow
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2019 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

arrows = findall(this.hParent, 'Tag', 'sliceomaticarrow');

pushset(arrows, 'FaceColor', [1, 0, 0]);
pushset(arrows, 'FaceAlpha', 0.2);

pushset(this.draggedArrow, 'FaceColor', [0, 1, 0]);
pushset(this.draggedArrow, 'FaceAlpha', 0.7);

hSlices = this.GetAllSlices();

for iSlice = 1:length(hSlices)
  fa = get(hSlices(iSlice), 'FaceAlpha');
  if isa(fa, 'double') && fa > 0.3
    pushset(hSlices(iSlice), 'FaceAlpha', 0.3);
    pushset(hSlices(iSlice), 'EdgeColor', 'none');
  else
    pushset(hSlices(iSlice), 'FaceAlpha', fa);
    pushset(hSlices(iSlice), 'EdgeColor', get(hSlices(iSlice), 'EdgeColor'));
  end
end

hIsoSurfs = this.GetAllIsos();

for iIso = 1:length(hIsoSurfs)
  fa = get(hIsoSurfs(iIso), 'FaceAlpha');
  if isa(fa, 'double') && fa > .3
    pushset(hIsoSurfs(iIso), 'FaceAlpha', .3);
    pushset(hIsoSurfs(iIso), 'EdgeColor', 'none');
  else
    pushset(hIsoSurfs(iIso), 'FaceAlpha', fa);
    pushset(hIsoSurfs(iIso), 'EdgeColor', get(hIsoSurfs(iIso), 'EdgeColor'));
  end
  hCap = getappdata(hIsoSurfs(iIso), 'isosurfacecap');
  if ~isempty(hCap)
    pushset(hCap, 'Visible', 'off');
  end
end

ss = getappdata(this.draggedArrow, 'arrowslice');

if isempty(ss)
  ss = getappdata(this.draggedArrow, 'arrowiso');
end

popset(ss, 'FaceAlpha');
popset(ss, 'EdgeColor');

wbmf = get(this.hFigure, 'WindowButtonMotionFcn');
% temporarily un-set figure mode
tempfiguremode(this.hFigure, 'unset');
% re-set callback deleted by figure-mode functions (work around Matlab bug)
set(this.hFigure, 'WindowButtonMotionFcn', wbmf);
pushset(this.hFigure, 'WindowButtonUpFcn', @(hObject, eventData) this.ArrowDragFinish());
pushset(this.hFigure, 'WindowButtonMotionFcn', @(hObject, eventData) this.callbacks('motion'));

% Doing this makes the tip invisible when visible is on.
showarrowtip(this.hParent, this.draggedArrow);

if ispc() && ~isMouseButtonPressed()
  % Call the window button up function now if the mouse button is no longer
  % pressed.
  this.ArrowDragFinish();
end

end
