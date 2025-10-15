function AxesDeleteFcn(this, hAxes, eventData)
%% Executes on deletion of main axes.
%
%       this.AxesDeleteFcn(hAxes, eventData, hParent)
%
% This function resets the menu and toolbar when the last sliceomatic is removed
% from the parent figure.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2018 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% figure might already be gone
if ~ishghandle(this.hParent)
  return
end

%% get structure with handles
d = getappdata(this.hParent, 'sliceomatic');

if isempty(d)
  return
end

%% remove objects
% remove sliceomatic from parent
rmappdata(this.hParent, 'sliceomatic');

% remove menus
delete_valid(this, 'fileMenu');
delete_valid(this, 'controlsMenu');
delete_valid(this, 'objectDefaultsMenu');
delete_valid(this, 'allSlicesMenu');
delete_valid(this, 'helpMenu');

% remove context menus
delete_valid(this, 'sliceContextMenu');
delete_valid(this, 'isoContextMenu');

% remove slider controls
delete_valid(this, 'hSliderX');
delete_valid(this, 'hSliderY');
delete_valid(this, 'hSliderZ');
delete_valid(this, 'hSliderIso');

% remove drop down boxes and titles
delete_valid(this, 'colormapPanel');
delete_valid(this, 'alphamapPanel');
delete_valid(this, 'orientationPanel');

% restore parent size changed function
if verLessThan('matlab', '8.4') % before R2014b
  if isfield(d, 'oldResizeFcn')
    set(this.hParent, 'ResizeFcn', d.oldResizeFcn);
  end
else
  if isfield(d, 'oldSizeChangedFcn')
    set(this.hParent, 'SizeChangedFcn', d.oldSizeChangedFcn);
  end
end

% For the last sliceomatic in figure, remove menubar and restore menubar and callbacks
if numel(findobj(this.hFigure, 'Tag', 'MainAxes')) <= 1
  % remove toolbars
  if exist('uitoolfactory', 'file') && isfield(d, 'toolbar') && ~isempty(d.toolbar) && ishghandle(d.toolbar)
    delete_valid(d, 'toolbar');
  end
  if isappdata(this.hFigure, 'sliceomatic_toolbar')
    rmappdata(this.hFigure, 'sliceomatic_toolbar');
  end
  try %#ok<TRYNC>
    cameratoolbar(this.hFigure, 'Hide');
  end
  if isappdata(this.hFigure, 'sliceomatic_oldtoolbar')
    set(this.hFigure, 'Toolbar', getappdata(d.handleFigure, 'sliceomatic_oldtoolbar'));
    rmappdata(this.hFigure, 'sliceomatic_oldtoolbar');
  end

  % reset callbacks
  if isappdata(this.hFigure, 'sliceomatic_oldbdf')
    oldbdf = getappdata(this.hFigure, 'sliceomatic_oldbdf');
    rmappdata(this.hFigure, 'sliceomatic_oldbdf');
    set(zoom(this.hFigure),     'ButtonDownFilter', oldbdf.zoom);
    set(pan(this.hFigure),      'ButtonDownFilter', oldbdf.pan);
    set(rotate3d(this.hFigure), 'ButtonDownFilter', oldbdf.rotate3d);
  end
  if isappdata(this.hFigure, 'sliceomatic_oldwbmf')
    set(this.hFigure, 'WindowButtonMotionFcn', getappdata(this.hFigure, 'sliceomatic_oldwbmf'));
    rmappdata(this.hFigure, 'sliceomatic_oldwbmf');
  end

  % un-register data cursor update function
  % FIXME: handle multiple sliceomatics in one figure
  DataCursorUpdateFcns = getappdata(this.hFigure, 'DataCursorUpdateFcns');
  if isfield(DataCursorUpdateFcns, 'sliceomatic')
    DataCursorUpdateFcns = rmfield(DataCursorUpdateFcns, 'sliceomatic');
    setappdata(this.hFigure, 'DataCursorUpdateFcns', DataCursorUpdateFcns);
  end
  set(this.hFigure, 'Pointer', 'arrow');

  % restore old menubar
  if isappdata(this.hFigure, 'sliceomatic_oldmenubar')
    set(this.hFigure, 'Menubar', getappdata(this.hFigure, 'sliceomatic_oldmenubar'));
  end

end

end

function delete_valid(this, fieldname)
%% delete object if it is valid

if ~isempty(this.(fieldname)) && ishghandle(this.(fieldname))
  delete(this.(fieldname));
end

end
