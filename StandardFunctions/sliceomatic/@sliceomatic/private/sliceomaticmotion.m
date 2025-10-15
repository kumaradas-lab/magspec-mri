function sliceomaticmotion(hFigure, eventData, hParent)
% Handle generic motion events for the figure window.

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc, 2015 Pure Devices

obj = hittest(hFigure);

% Some objects get a special pointer when the mouse waves over them.
% Get it from appdata.
if ~isempty(obj)
  t = getappdata(obj, 'motionpointer');
  cc = get(hFigure, 'Pointer');

  if t
    newc = t;
  else
    newc = get(0, 'DefaultFigurePointer');
  end

  if isa(newc, 'char') && isa(cc, 'char') && ~strcmp(newc, cc)
    setpointer(hFigure, newc);
  end
end

if ~ishghandle(hParent)
  return;
end

d = getappdata(hParent, 'sliceomatic');

if ~isfield(d, 'displayComplete') || ~d.displayComplete
  % window is not drawn yet. do nothing
  return
end

% The motion meta slice is a line that is managed here.
% There is only one.
if isempty(d.motionmetaslice)
  d.motionmetaslice(1) = line('Parent', d.axmain, ...
    'Visible', 'off', ...
    'LineStyle', '-', ...
    'Marker', 'none', ...
    'LineWidth', 2, ...
    'Clipping', 'off', ...
    'Color', 'w');
  d.motionmetaslice(2) = line('Parent', d.axmain, ...
    'Visible', 'off', ...
    'LineStyle', '--', ...
    'Marker', 'none', ...
    'LineWidth', 2, ...
    'Clipping', 'off', ...
    'Color', 'k');
  setappdata(hParent, 'sliceomatic', d);
end

showarrowtip(hParent, obj);

if isempty(obj) || (obj ~= d.axx && obj ~= d.axy && obj ~= d.axz)
  if ishandle(d.motionmetaslice)
    set(d.motionmetaslice, 'Visible', 'off');
  end
  return
end

% OBJ can only be an Axes because of the previous IF statement.
aa = obj;
apos = get(aa, 'CurrentPoint');

xl = d.xlim;
yl = d.ylim;
zl = d.zlim;

if aa==d.axx || aa==d.axiso
  if aa==d.axiso
    % eh?
    return
  else
    xdata = [ apos(1,1) apos(1,1) apos(1,1) apos(1,1) apos(1,1) ];
    ydata = [ yl(1) yl(2) yl(2) yl(1) yl(1) ];
    zdata = [ zl(2) zl(2) zl(1) zl(1) zl(2) ];
  end
else
  % We are moving a Y or Z slice
  if aa==d.axy
    ydata = [ apos(1,2) apos(1,2) apos(1,2) apos(1,2) apos(1,2) ];
    xdata = [ xl(1) xl(2) xl(2) xl(1) xl(1) ];
    zdata = [ zl(2) zl(2) zl(1) zl(1) zl(2) ];
  else
    zdata = [ apos(1,2) apos(1,2) apos(1,2) apos(1,2) apos(1,2) ];
    ydata = [ yl(1) yl(2) yl(2) yl(1) yl(1) ];
    xdata = [ xl(2) xl(2) xl(1) xl(1) xl(2) ];
  end
end

%set(d.motionmetaslice(1), 'Visible', 'on', ...
%  'XData', xdata, 'YData', ydata, 'ZData', zdata);
set(d.motionmetaslice, 'Visible', 'on', ...
  'XData', xdata, 'YData', ydata, 'ZData', zdata);

end
