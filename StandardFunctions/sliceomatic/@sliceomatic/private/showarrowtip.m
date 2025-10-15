function showarrowtip(hParent, arrow)
%% Display a small text field next to the ARROW.
% Depends on tipdata being set on the handle to ARROW.
%
% Copyright 2000, 2001, 2002, 2003, 2004, 2005 The MathWorks Inc
%
% ------------------------------------------------------------------------------
% (C) Copyright 2019-2021 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

%%
d = getappdata(hParent, 'sliceomatic');

if isempty(arrow) || ~ishghandle(arrow)
  return;
end

ctrlarrow = getappdata(arrow, 'controlarrow');

if ~isempty(ctrlarrow) && ishghandle(ctrlarrow)
  % In this case, the slice or isosurface passed in most likely from the motion
  % callback.  Lets redirect our input argument and show the tip anyway.
  
  % See the "arrow" private function for the setting of this value;
  arrow = ctrlarrow;
end

tipdata = getappdata(arrow, 'tipdata');

if ~isempty(tipdata) && ~isemptyfield(tipdata, 'value')
  if abs(tipdata.value) < 1e-2 && tipdata.value ~= 0
    tipDataString = sprintf('Value: %1.3e', tipdata.value);
  else
    tipDataString = sprintf('Value: %1.3f', tipdata.value);
  end
  set(d.tip, ...
    'Parent', tipdata.parentaxes, ...
    'String', tipDataString, ...
    'Units', 'data', ...
    'Position', tipdata.position, ...
    'VerticalAlignment', tipdata.verticalalign, ...
    'HorizontalAlignment', tipdata.horizontalalign);
  set(d.tip, 'Units', 'pixels');
  set(d.tip, 'Visible', 'on');
elseif ishandle(d.tip)
  set(d.tip, 'Visible', 'off');
end

end
