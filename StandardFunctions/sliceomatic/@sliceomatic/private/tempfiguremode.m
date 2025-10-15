function tempfiguremode(hFigure, action)
%% Un-set or re-set the figure mode
%
%       tempfiguremode(hFigure, action)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2019 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

persistent currentMode

switch action
  case 'unset'
    hManager = uigetmodemanager(hFigure);
    currentMode = hManager.CurrentMode;
    zoom(hFigure, 'off');
    pan(hFigure, 'off');
    rotate3d(hFigure, 'off');
    datacursormode(hFigure, 'off');
  case 'reset'
    if isempty(currentMode)
      return
    end
    switch currentMode.Name
      case 'Exploration.Zoom'
        zoom(hFigure, 'on');
      case 'Exploration.Pan'
        pan(hFigure, 'on');
      case 'Exploration.Rotate3d'
        rotate3d(hFigure, 'on');
    end
end

end
