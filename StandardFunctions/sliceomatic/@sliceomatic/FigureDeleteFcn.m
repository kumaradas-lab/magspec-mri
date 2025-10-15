function FigureDeleteFcn(this, hFigure, eventData)
%% Executes on deletion of parent figure.
%
%       this.FigureDeleteFcn(hFigure, eventData)
%
% The parent figure might be deleted after sliceomatic was re-parented.
% In this case, all objects that are bound to the figure have to be re-parented
% before they become invalid.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% Check if new figure is valid
newFigure = ancestor(this.hAxes, 'figure');

if isempty(newFigure) || ~ishghandle(newFigure, 'figure')
  return;
end

this.hFigure = newFigure;

%% Figure menus
set(this.fileMenu, 'Parent', this.hFigure);
set(this.controlsMenu, 'Parent', this.hFigure);
set(this.objectDefaultsMenu, 'Parent', this.hFigure);
set(this.allSlicesMenu, 'Parent', this.hFigure);
set(this.helpMenu, 'Parent', this.hFigure);

%% Context menus
set(this.sliceContextMenu, 'Parent', this.hFigure);
set(this.isoContextMenu, 'Parent', this.hFigure);

%% Set this delete function to the new figure.
set(this.hFigure, 'DeleteFcn', @this.FigureDeleteFcn);  % FIXME: Should we care to restore a previous delete function when we are done?

%% Set up WindowButtonMotionFcn
if ~isappdata(this.hFigure, 'sliceomatic_oldwbmf')
  oldwbmf = get(this.hFigure, 'WindowButtonMotionFcn');
  setappdata(this.hFigure, 'sliceomatic_oldwbmf', oldwbmf);
  if ~isempty(oldwbmf) && ~(isequal(oldwbmf, @sliceomaticmotion) || ...
      (iscell(oldwbmf) && isequal(oldwbmf{1}, @sliceomaticmotion)))
    warning('sliceomatic:OverwriteWindowButtonMotionFcn', ...
      'WindowButtonMotionFcn is overwritten by sliceomatic');
  end
  set(this.hFigure, 'WindowButtonMotionFcn', {@sliceomaticmotion, this.hParent});
end

%% Button Down Filter Functions
set(zoom(this.hFigure),     'ButtonDownFilter', @(obj,eventdata) ButtonDownFilterFcn(this.hParent));
set(pan(this.hFigure),      'ButtonDownFilter', @(obj,eventdata) ButtonDownFilterFcn(this.hParent));
set(rotate3d(this.hFigure), 'ButtonDownFilter', @(obj,eventdata) ButtonDownFilterFcn(this.hParent));

end
