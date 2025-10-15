function isocontrols(fig, onoff)
% Set up FIG to have an ISO surface controller on the bottom.
% ONOFF indicates if the controller is being turned ON or OFF

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% Check variables
narginchk(2, 2);

d = getappdata(fig, 'sliceomatic');

if onoff
  lim = [min(d.data(isfinite(d.data))), max(d.data(isfinite(d.data)))];
  if isempty(lim)
    lim = [-1, 1];
  end
  if lim(1) == lim(2)
    lim(1) = 0.9*lim(1);
    lim(2) = 1.1*lim(2);
  end

  set(d.axiso, 'HandleVisibility', 'on');
  set(d.axiso, ...
    'XLim', lim, ...
    'YLim', [1, 5], ...
    'CLim', [1, 64]);
  l = activelabel(d.axiso, 'Title', 'Iso Surface Controller');

  pos_axis = get(d.axiso, 'OuterPosition');
  set(l, 'Units', 'normalized');
  pos_label = get(l, 'Position');
  pos_label(1) = pos_axis(1) + pos_axis(3)/2;
  pos_label(2) = pos_axis(2) - .6;
  set(l, 'Position', pos_label, ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'center');

  image('Parent', d.axiso, 'CData', 1:64, 'CDataMapping', 'scaled', ...
    'XData', lim, 'YData', [1, 5], 'Tag', 'sliceomaticisocontrolimage', ...
    'HitTest', 'off');
  set(d.axiso, 'HandleVisibility', 'off');

else
  % Turn off the controller
  delete(findobj(d.axiso, 'Type', 'image'));
end

end
