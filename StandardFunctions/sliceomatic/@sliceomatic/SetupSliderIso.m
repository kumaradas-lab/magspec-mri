function SetupSliderIso(this, onoff)
%% Add a slider for iso surfaces to the bottom of the figure
% ONOFF indicates if the controller is being turned ON or OFF

% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

d = getappdata(this.hParent, 'sliceomatic');

if onoff
  lim = [min(d.data(isfinite(d.data))), max(d.data(isfinite(d.data)))];
  if isempty(lim) || all(lim == 0) || all(isnan(lim))
    lim = [-1, 1];
  end
  if lim(1) == lim(2)
    if lim(1) > 0
      lim(1) = 0.9*lim(1);
      lim(2) = 1.1*lim(2);
    else
      lim(1) = 1.1*lim(1);
      lim(2) = 0.9*lim(2);
    end
  end

  set(this.hSliderIso, 'HandleVisibility', 'on');
  set(this.hSliderIso, ...
    'XLim', lim, ...
    'YLim', [1, 5], ...
    'CLim', [1, 64]);
  l = activelabel(this.hSliderIso, 'Title', 'Iso Surface Controller');

  pos_axis = get(this.hSliderIso, 'OuterPosition');
  set(l, 'Units', 'normalized');
  pos_label = get(l, 'Position');
  pos_label(1) = pos_axis(1) + pos_axis(3)/2;
  pos_label(2) = pos_axis(2) - .6;
  set(l, 'Position', pos_label, ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'center');

  this.hSliderIsoImage = image(...
    'Parent', this.hSliderIso, 'CData', 1:64, 'CDataMapping', 'scaled', ...
    'XData', lim, 'YData', [1, 5], 'Tag', 'sliceomaticisocontrolimage', ...
    'HitTest', 'off');
  set(this.hSliderIso, 'HandleVisibility', 'off');

else
  % Turn off the controller
  delete(findobj(this.hSliderIso, 'Type', 'image'));
end

end
