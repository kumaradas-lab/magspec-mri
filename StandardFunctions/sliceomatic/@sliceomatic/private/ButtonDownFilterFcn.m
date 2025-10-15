function modeOverwrite = ButtonDownFilterFcn(handleParent)
%% Filter function for figure modes
%
% modeOverwrite [output] a logical flag to determine whether the
%               pan/zoom/rotate operation should take place (for
%               'modeOverwrite' set to 'false') 
%               or the 'ButtonDownFcn' property of the object should 
%               take precedence (when 'modeOverwrite' is 'true')
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if ~ishghandle(handleParent)
  modeOverwrite = false;
  return;
end

hFigure = ancestor(handleParent, 'figure');
obj = hittest(hFigure);
% find all sliceomatics in current figure
handleParents = get(findobj(hFigure, 'Tag', 'MainAxes'), 'Parent');
numHandleParents = numel(handleParents);
if numHandleParents > 1
  for iSlix = numel(handleParents):-1:1
    d(iSlix) = getappdata(handleParents{iSlix}, 'sliceomatic');
  end
else
  d = getappdata(handleParents, 'sliceomatic');
end
controlArrowHandle = 0;
try
    [controlArrowHandle, surfHandle] = getarrowslice();
    if isempty(surfHandle)
        % getarrowslice returns gco for any object
        controlArrowHandle = 0;
    end
catch ME
    disp(getReport(ME))
end

modeOverwrite = false;
if ~isempty(obj) && ...
        (any(obj == [d(:).axx]) || any(obj == [d(:).axy]) || ...
         any(obj == [d(:).axz]) || any(obj == [d(:).axiso]) || ...
        any(obj == controlArrowHandle))
    % overwrite pan, zoom or rotate3d mode when one of the sliders or
    % arrows are clicked
    modeOverwrite = true;
end

end
