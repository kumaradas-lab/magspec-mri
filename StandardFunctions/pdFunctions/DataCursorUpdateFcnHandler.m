function outputTxt = DataCursorUpdateFcnHandler(pointDataTip, eventObj)
%% Handler for figures that provide multiple data cursor update functions
%
% To register, set a (unique) field with a function handle in the structure
% stored as the 'DataCursorUpdateFcns' appdata of the figure. The function
% handle must be callable with [outputTxt] = @(pointDataTip, eventObj).
% The first function that return a non-empty 
%
% pointDataTip A PointDataTip object (undocumented Matlab and currently not used)
% eventObj     Object containing event data structure
% outputTxt    Data cursor text (string or cell array of strings)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% call all DataCursorUpdateFcns until we get a valid string
hFigure = ancestor(get(eventObj, 'Target'), 'figure');

DataCursorUpdateFcns = getappdata(hFigure, 'DataCursorUpdateFcns');

% default if no data cursor update functions are registered (any more)
outputTxt = '';
if isstruct(DataCursorUpdateFcns)
  fields = fieldnames(DataCursorUpdateFcns);
  for iField = 1:numel(fields)
    outputTxt = DataCursorUpdateFcns.(fields{iField})(pointDataTip, eventObj);
    if ~isempty(outputTxt), break; end
  end
end

if isempty(outputTxt)
  % In case some axes didn't "register" a dedicated data cursor update function,
  % provide some default behavior.
  pos = get(eventObj, 'Position');
  if numel(pos) < 3
    outputTxt = {['X: ' num2str(pos(1), 3)], ['Y: ' num2str(pos(2), 3)]};
  else
    outputTxt = {['X: ' num2str(pos(1), 3)], ['Y: ' num2str(pos(2), 3)], ['Z: ' num2str(pos(3), 3)]};
  end
  target = get(eventObj, 'Target');
  switch get(target, 'Type')
    case 'image'
      cdata = get(target, 'CData');
      xdata = get(target, 'XData');
      if numel(xdata) == 2
        idx_x = pos(1);
      else
        diff_x = abs(xdata-pos(1));
        idx_x = find(diff_x == min(diff_x), 1);
      end
      ydata = get(target, 'YData');
      if numel(ydata) == 2
        idx_y = pos(2);
      else
        diff_y = abs(ydata-pos(2));
        idx_y = find(diff_y == min(diff_y), 1);
      end
      outputTxt{end+1} = ['V: ' num2str(cdata(idx_y, idx_x), 3)];
  end
end

end

