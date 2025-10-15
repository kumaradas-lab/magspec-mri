function d = figtoolbar(d)
% Set up the toolbar for Sliceomatic within structure D

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc

if ~isappdata(d.handleFigure, 'sliceomatic_oldtoolbar')
 setappdata(d.handleFigure, 'sliceomatic_oldtoolbar', get(d.handleFigure, 'Toolbar'));
end
set(d.handleFigure, 'Toolbar', 'none');

if exist('uitoolfactory', 'file')
  % Create a toolbar with just the elements useful on sliceomatic
  % on it.

  if isappdata(d.handleFigure, 'sliceomatic_toolbar')
    d.toolbar = getappdata(d.handleFigure, 'sliceomatic_toolbar');
  else
    d.toolbar = uitoolbar('Parent', d.handleFigure);
    uitoolfactory(d.toolbar, 'Annotation.InsertRectangle');
    uitoolfactory(d.toolbar, 'Annotation.InsertEllipse');
    uitoolfactory(d.toolbar, 'Annotation.InsertTextbox');
    uitoolfactory(d.toolbar, 'Annotation.InsertArrow');
    uitoolfactory(d.toolbar, 'Annotation.InsertLine');
    uitoolfactory(d.toolbar, 'Exploration.ZoomIn');
    uitoolfactory(d.toolbar, 'Exploration.ZoomOut');
    uitoolfactory(d.toolbar, 'Exploration.Pan');
    uitoolfactory(d.toolbar, 'Exploration.Rotate');
    uitoolfactory(d.toolbar, 'Exploration.DataCursor');
    setappdata(d.handleFigure, 'sliceomatic_toolbar', d.toolbar);
  end

  cameratoolbar(d.handleFigure, 'Show');
  cameratoolbar(d.handleFigure, 'ToggleSceneLight');
else
  % We are in R13 or earlier
  try
    cameratoolbar(d.handleFigure, 'Show');
    cameratoolbar(d.handleFigure, 'ToggleSceneLight');
    %cameratoolbar('SetMode','orbit');
  catch
    disp('Could not display the camera toolbar.');
  end
end

end
