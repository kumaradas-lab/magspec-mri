function [hAxes, config] = plot_data_1D(HW, data_1D, hAxes, raiseWindow, config)
%% Plot RF signal (abs, real, imag), its phase and its frequency offset to AQ frequency
%
%       [hAxes, config] = plot_data_1D(HW, data_1D, hAxes, raiseWindow, config)
%
% INPUT:
%   HW            HW structure or object
%   data_1D       structure with measurement data
%   hAxes         array with graphics handles for the three axes containing the
%                 RF signal (abs, real, imag), its phase and its frequency
%                 offset to AQ frequency (see corresponding output argument).
%                 Or: graphics handle to a valid parent for these three axes
%                 (uipanel or figure handle).
%                 If omitted, empty or 1, figure 100 is used.
%   raiseWindow   Boolean value. If false and hAxes is not empty, the figure
%                 window does not steal the focus if it already exists.
%                 (Default: true)
%   config        An optional structure with the following fields. If the these
%                 fields are omitted or empty, default values are used:
%     hParent       Handle to a figure or uipanel that is used as a parent for
%                   the axes (default: 100).
%     figureTitle   If "hParent" is a figure, this string is set as the figure
%                   name (default: 'Timeline continuous').
%     timeFieldname The fieldname in the structure "data_1D" that is used as the
%                   timeline. Valid values are "time_all" and "time_of_tRep"
%                   (default: "time_all").
%     plotData      Boolean value. If true, plot the absolute value (and its
%                   real and imaginary parts, see "plotDataRealImag" below).
%                   (Default: true)
%     plotDataRealImag
%                   Boolean value. If true, the real and imaginary parts of the
%                   signal is plotted additionally to its absolute value (see
%                   "plotData" above). (Default: true)
%     plotPhase     Boolean value. If true, plot the phase of the data (default:
%                   true).
%     plotFreq      Boolean value. If true, plot the frequency (offset) of the
%                   data (default: true).
%     plotDataYLim  2-element vector with the y-limits in Tesla for the data
%                   plot. If this is "auto", the y-limits mode is set to "auto"
%                   instead (default).
%     omitnan       Boolean value. If true, connect points in plot even if they
%                   are separated by NaN values (default: false).
%
% OUTPUT:
%   hAxes         array with graphics handles for the three axes containing the
%                 RF signal (abs, real, imag), its phase and its frequency
%                 offset to AQ frequency (see corresponding input argument).
%   config        Same as the input argument with the actually used settings.
%
% ------------------------------------------------------------------------
% (C) Copyright 2011-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------

%% Default input
if nargin < 3, hAxes = []; end
if nargin < 4, raiseWindow = true; end
if nargin < 5, config = []; end

config = set_EmptyField(config, 'hParent', 100);
config = set_EmptyField(config, 'figureTitle', 'Timeline continuous');
config = set_EmptyField(config, 'timeFieldname', 'time_all');
config = set_EmptyField(config, 'plotData', true);
config = set_EmptyField(config, 'plotPhase', true);
config = set_EmptyField(config, 'plotFreq', true);
config = set_EmptyField(config, 'plotDataRealImag', true);
config = set_EmptyField(config, 'iData', 1:numel(data_1D));
if ~isfield(config, 'plotDataYLim'), config.plotDataYLim = []; end
if ~isnumeric(config.plotDataYLim) || numel(config.plotDataYLim) ~= 2
  config.plotDataYLim = [];
end
config = set_EmptyField(config, 'omitnan', false);

%% handle multiple data channels
if ~isscalar(data_1D)
  configIn = config;
  hAxesIn = hAxes;
  for iData = flip(configIn.iData(:).')  % master last so it is on top
    if iscell(hAxesIn), hAxes = hAxesIn{iData}; else hAxes = hAxesIn; end
    if numel(configIn.hParent) >= iData
      config.hParent = configIn.hParent(iData);
    else
      if ishghandle(configIn.hParent(1), 'figure') || isnumeric(configIn.hParent(1))
        config.hParent = double(configIn.hParent) + iData-1;
      else
        config.hParent = 100 + iData-1;
      end
    end
    config.figureTitle = sprintf('%s (data channel #%d)', configIn.figureTitle, iData);
    config.iData = iData;
    [hAxesOut{iData}, configOut(iData)] = plot_data_1D(HW, data_1D(iData), hAxes, raiseWindow, config);
  end
  hAxes =  hAxesOut;
  config = configOut;
  return;
end

%% "multi-threading protection"
% block execution while function is running
function oldState = isRunning(newState)
  persistent running
  if isempty(running), running = false; end
  if nargout > 0, oldState = running; end
  if nargin > 0, running = newState; end
end

if isRunning(), return; end
isRunning(true);
cleanupObj = onCleanup(@() isRunning(false));

%% handle input argument "hAxes"
plotAxes = logical([config.plotData, config.plotPhase, config.plotFreq]);
if ~any(plotAxes)
  return;
end
lastAxes = find(plotAxes, 1, 'last');

% check type of hAxes and set graphic handles accordingly
if numel(hAxes) == 1 && ...
    (ishghandle(hAxes, 'uipanel') || (ishghandle(hAxes, 'figure') || ...
     (isa(hAxes, 'double') && mod(hAxes, 1) == 0)) && hAxes ~= 1)
  hParent = hAxes;
  hAxes = [];
elseif ~isempty(hAxes) && numel(hAxes) >= lastAxes && ...
    all(ishghandle(hAxes(plotAxes(1:lastAxes)), 'axes'))
  hParent = get(hAxes(find(ishghandle(hAxes, 'axes'), 1, 'first')), 'Parent');
else
  hParent = config.hParent;
end

hFigure = ancestor(hParent, 'figure');
if isempty(hFigure), hFigure = hParent; end % hParent is handle to a non-existent figure

% Is the (future?) parent the figure directly?
if ishghandle(hParent, 'figure') || (isa(hParent, 'double') && mod(hParent, 1) == 0)
  isFigure = true;
else
  isFigure = false;
end
if isFigure && ~ishghandle(hFigure, 'figure')
  % create figure if not yet present
  hFigure = figure(hFigure);
end

%% find existing axes
if ishghandle(hParent) && (isempty(hAxes) || ~all(ishghandle(hAxes, 'axes')))
  hAxes = findobj(get(hParent, 'Children'), 'flat', '-regexp', 'Tag', 'plotData1D_Axes[1-3]');
end

% sort axes according to "Tag"s
hAxesTmp = hAxes(ishghandle(hAxes, 'axes'));
hAxes = [];
axesTag = get(hAxesTmp, 'Tag');
if numel(hAxesTmp) > 1
  idx = cellfun(@(str) str2double(str(end)), axesTag);
  hAxes(idx) = hAxesTmp;
elseif numel(hAxesTmp) == 1
  hAxes(str2double(axesTag(end))) = hAxesTmp;
end

%% setup data cursor update function
if isempty(hAxes)
  % assume that this is the first time we use this figure
  DataCursorUpdateFcns = getappdata(hFigure, 'DataCursorUpdateFcns');
  % FIXME: The last argument to the data cursor function should be the device index
  iDevice = data_1D.device;
  DataCursorUpdateFcns.plotData1D = @(pointDataTip, eventData) ...
    plotData1D_DataCursorFcn(pointDataTip, eventData, HW, iDevice);  
  setappdata(hFigure, 'DataCursorUpdateFcns', DataCursorUpdateFcns);
  hdcm = datacursormode(hFigure);
  set(hdcm, 'UpdateFcn', @DataCursorUpdateFcnHandler);

  set(hParent, 'DeleteFcn', @plotData1D_DeleteFcn);
end

%% clear foreign elements from parent
hKids = get(hParent, 'Children');
% FIXME: Annotation text box sometimes (when?) switches 'Tag'
hForeign = findobj(hKids, 'flat', '-not', {'-regexp', 'Tag', 'plotData1D_Axes[1-3]', ...
                                           '-or', 'Tag', 'plotData1D_legend', ...
                                           '-or', 'Type', 'annotationpane'});
delete(hForeign);

%% Check if hAxes is (still) valid and re-create axes if necessary
if isempty(hAxes) || numel(hAxes) < lastAxes || any(~ishghandle(hAxes(plotAxes(1:lastAxes)), 'axes'))
  fig_title = config.figureTitle;
  found_axes(1:3) = false;
  if ~raiseWindow && ishghandle(hParent, 'figure')
    % find axes in parent (= figure)
    hKids = get(hParent, 'Children');
    for i_axes = find(plotAxes)
      hAxis = findobj(hKids, 'flat', 'Tag', sprintf('plotData1D_Axes%d', i_axes));
      if ishghandle(hAxis, 'axes')
        hAxes(i_axes) = hAxis;
        found_axes(i_axes) = true;
      end
    end
    if ~all(found_axes | ~plotAxes)
      % some axes are missing, delete all existing children
      % FIXME: Re-use existing axes. Only re-create missing ones.
      if isFigure
        hParent = clf(hParent);
      else
        delete(get(hParent, 'Children'));
      end
    end
  else
    % raise figure and/or delete all existing children
    if ishghandle(hParent, 'figure')
      % legacy behavior
      hParent = clf(hParent);
    elseif ishghandle(hParent)
      % probably inside a panel
      delete(get(hParent, 'Children'));
    else
      hParent = figure(hParent);
    end
    hAxes(:) = 0;
  end
  if isFigure
    set(hFigure, 'Name', fig_title, 'NumberTitle', 'off');
  end
  if ~all(found_axes | ~plotAxes)
    nAxes = 0;
    if plotAxes(1)
      nAxes = nAxes + 1;
      hAxes(1) = subplot(sum(plotAxes), 1, nAxes, 'Parent', hParent, 'Tag', 'plotData1D_Axes1');
    end
    if plotAxes(2)
      nAxes = nAxes + 1;
      hAxes(2) = subplot(sum(plotAxes), 1, nAxes, 'Parent', hParent, 'Tag', 'plotData1D_Axes2');
    end
    if plotAxes(3)
      nAxes = nAxes + 1;
      hAxes(3) = subplot(sum(plotAxes), 1, nAxes, 'Parent', hParent, 'Tag', 'plotData1D_Axes3');
    end
    if numel(hAxes(plotAxes)) > 1
      linkaxes(hAxes(plotAxes), 'x');
    end
    if isFigure
      set(zoom(hFigure), 'Motion', 'both', 'Enable', 'on');
    end
  end
end

%% clear axes that are no longer needed
numAxesChanged = false;
if numel(hAxes) > lastAxes || any(ishghandle(hAxes(~plotAxes(1:numel(hAxes))), 'axes'))
  numAxesChanged = true;
  hDelAxes = hAxes(~plotAxes(1:numel(hAxes)));
  delete(hDelAxes(ishghandle(hDelAxes, 'axes')));
  hAxes(~plotAxes(1:numel(hAxes))) = 0;
  % re-position remaining axes
  for hCurrentAxes = hAxes(plotAxes)
    if ishghandle(hCurrentAxes, 'axes')
      subplot(sum(plotAxes), 1, sum(plotAxes(1:find(hCurrentAxes == hAxes))), hCurrentAxes);
    end
  end
end

hParent = get(hAxes(find(ishghandle(hAxes, 'axes'), 1, 'first')), 'Parent');

%% set focus to figure
hFigure = ancestor(hAxes(find(ishghandle(hAxes, 'axes'), 1, 'first')), 'figure');
if raiseWindow
  figure(hFigure);
end

%% plot or update figure
persistent hLabelY hMeanFreq
isNewAxes = false;
%% data
if config.plotData && ~isempty(data_1D.data)
  ylabelStr = {HW.RX(data_1D.device).AmplitudeName; ['in ', HW.RX(data_1D.device).AmplitudeUnit]};
  singleLine = all(imag(data_1D.data)==0) || ~config.plotDataRealImag;

  % search graphics objects
  hKids = get(hAxes(1), 'Children');
  found_lines(1:3) = false;
  h_lines = zeros(1,3);
  for i_line = 1:3
    h_line = findobj(hKids, 'flat', 'Tag', sprintf('plotData1D_Line%d', i_line));
    if ishghandle(h_line, 'line')
      h_lines(i_line) = h_line;
      found_lines(i_line) = true;
    end
  end

  if config.omitnan
    plotIdx = ~isnan(data_1D.(config.timeFieldname)) & ~isnan(data_1D.data);
    % But don't omit NaNs at wraps
    t = diff(plotIdx);
    lastValids = find([t;0] == -1);
    firstValids = find([0;t] == 1);
    if firstValids(1) < lastValids(1), firstValids(1) = false; end
    lastValids = lastValids(1:numel(firstValids));
    isWrap = data_1D.(config.timeFieldname)(firstValids) < ...
      data_1D.(config.timeFieldname)(lastValids);
    idxWrap = arrayfun(@(x,y) x:y, lastValids(isWrap), firstValids(isWrap), 'UniformOutput', false);
    idxWrap = horzcat(idxWrap{:});
    plotIdx(idxWrap) = true;
  else
    plotIdx = 1:numel(data_1D.(config.timeFieldname));
  end
  if ~all(found_lines)
    if singleLine
      nanData = nan(size(data_1D.data(plotIdx)));
      if ~config.plotDataRealImag
        h_lines(1) = plot(hAxes(1), ...
          data_1D.(config.timeFieldname)(plotIdx), ...
          abs(data_1D.data(plotIdx))/HW.RX(data_1D.device).AmplitudeUnitScale, ...
          'Tag', 'plotData1D_Line1');
        hold(hAxes(1), 'all');
        h_lines(2) = plot(hAxes(1), ...
          data_1D.(config.timeFieldname)(plotIdx), nanData, ...
          'Tag', 'plotData1D_Line2');
      else
        h_lines(1) = plot(hAxes(1), ...
          data_1D.(config.timeFieldname)(plotIdx), nanData, ...
          'Tag', 'plotData1D_Line1');
        hold(hAxes(1), 'all');
        h_lines(2) = plot(hAxes(1), ...
          data_1D.(config.timeFieldname)(plotIdx), ...
          real(data_1D.data(plotIdx))/HW.RX(data_1D.device).AmplitudeUnitScale, ...
          'Tag', 'plotData1D_Line2');
      end
      h_lines(3) = plot(hAxes(1), ...
        data_1D.(config.timeFieldname)(plotIdx), nanData, ...
        'Tag', 'plotData1D_Line3');
      hold(hAxes(1), 'off');
    else
      h_lines(1) = plot(hAxes(1), ...
        data_1D.(config.timeFieldname)(plotIdx), ...
        abs(data_1D.data(plotIdx))/HW.RX(data_1D.device).AmplitudeUnitScale, ...
        'Tag', 'plotData1D_Line1');
      hold(hAxes(1), 'all');
      h_lines(2) = plot(hAxes(1), ...
        data_1D.(config.timeFieldname)(plotIdx), ...
        real(data_1D.data(plotIdx))/HW.RX(data_1D.device).AmplitudeUnitScale, ...
        'Tag', 'plotData1D_Line2');
      h_lines(3) = plot(hAxes(1), ...
        data_1D.(config.timeFieldname)(plotIdx), ...
        imag(data_1D.data(plotIdx))/HW.RX(data_1D.device).AmplitudeUnitScale , ...
        'Tag', 'plotData1D_Line3');
      hold(hAxes(1), 'off');
    end
    numAQs = sum(isnan(data_1D.(config.timeFieldname)));
    if length(data_1D.(config.timeFieldname))/numAQs == 2 && ~config.omitnan
      set(h_lines, 'LineStyle', 'none', 'Marker', '.');
    else
      set(h_lines, 'LineStyle', '-', 'Marker', 'none');
    end
    set(hAxes(1), 'Tag', 'plotData1D_Axes1');
    lh = legend(hAxes(1), {'abs', 'real', 'imag'}, 'Orientation', 'horizontal', ...
      'Location', 'north', 'Tag', 'plotData1D_legend');
    legend(hAxes(1), 'boxoff');
    % move legend outside manually to keep the axes height
    legendPos = get(lh, 'Position');
    axPos = get(hAxes(1), 'Position');
    legendPos(2) = axPos(2) + axPos(4);
    set(lh, 'Position', legendPos);
    if singleLine
      set(lh, 'Visible', 'off');
    end
    if ~verLessThan('Matlab', '9.2')
      set(lh, 'AutoUpdate', 'off'); % Do not automatically add new entries to legend
    end
    grid(hAxes(1), 'on');
    isNewAxes = true;
    ylabel(hAxes(1), ylabelStr);
    if isempty(config.plotDataYLim)
      set(hAxes(1), 'YLimMode', 'auto');
    else
      set(hAxes(1), 'YLim', config.plotDataYLim/HW.RX(data_1D.device).AmplitudeUnitScale);
    end
  else
    lh = findobj(hParent, 'Tag', 'plotData1D_legend');
    if singleLine
      nanData = nan(size(data_1D.data(plotIdx)));
      if ~config.plotDataRealImag
        yData = {abs(data_1D.data(plotIdx))/HW.RX(data_1D.device).AmplitudeUnitScale; nanData; nanData};
      else
        yData = {nanData; real(data_1D.data(plotIdx))/HW.RX(data_1D.device).AmplitudeUnitScale; nanData};
      end
      set(lh, 'Visible', 'off');
    else
      yData = mat2cell([abs(data_1D.data(plotIdx)), real(data_1D.data(plotIdx)), imag(data_1D.data(plotIdx))] ./ ...
                       HW.RX(data_1D.device).AmplitudeUnitScale, size(data_1D.data(plotIdx),1), [1 1 1]).';
      set(lh, 'Visible', 'on');
    end
    set(h_lines, 'XData', data_1D.(config.timeFieldname)(plotIdx), ...
      {'YData'}, yData);
    numAQs = sum(isnan(data_1D.(config.timeFieldname)));
    if length(data_1D.(config.timeFieldname))/numAQs==2 && ~config.omitnan
      set(h_lines, 'LineStyle', 'none', 'Marker', '.');
    else
      set(h_lines, 'LineStyle', '-', 'Marker', 'none');
    end

    if ~isequal(get(hAxes(1), 'YLabel'), ylabelStr)
      ylabel(hAxes(1), ylabelStr);
    end
  end
end

%% phase
if config.plotPhase && ~isempty(data_1D.data)
  % cla(haxes(2),'reset')
  h_lines = get(hAxes(2), 'Children');
  if numel(h_lines) > 0
    h_lines(~strcmpi(get(h_lines, 'Type'), 'line')) = [];
  end
  if numel(h_lines) ~= 1
    h_lines = plot(hAxes(2), data_1D.(config.timeFieldname), angle(data_1D.data), 'Tag', 'plotData1D_LinePhase');
    set(hAxes(2), 'Tag', 'plotData1D_Axes2');
    grid(hAxes(2), 'on');
    isNewAxes = true;
    ylabel(hAxes(2), {'Phase'; 'in rad'});
  else
    set(h_lines, 'XData', data_1D.(config.timeFieldname), ...
      'YData', angle(data_1D.data));
  end
  if length(data_1D.(config.timeFieldname))/sum(isnan(data_1D.(config.timeFieldname)))==2 && ~config.omitnan
    set(h_lines, 'LineStyle', 'none', 'Marker', '.');
  else
    set(h_lines, 'LineStyle', '-', 'Marker', 'none');
  end
end

%% frequency (offset)
if config.plotFreq && ~isempty(data_1D.data)
  % cla(haxes(3), 'reset')
  h_lines = get(hAxes(3), 'Children');
  if numel(h_lines) > 0
    h_lines(~strcmpi(get(h_lines, 'Type'), 'line')) = [];
  end
  grad_time = gradient(data_1D.(config.timeFieldname));
  offsetFrequency = gradient(unwrap(angle(data_1D.data)))./2./pi./grad_time;
  signalFrequency = data_1D.AqFrequency - offsetFrequency;
  iv = ~isnan(data_1D.data);
  meanOffset = -meanNAN(offsetFrequency(iv).*abs(data_1D.data(iv))/sum(abs(data_1D.data(iv))))*sum(iv);
  if any(data_1D.AqFrequency(~isnan(data_1D.AqFrequency)) ~= ...
      data_1D.AqFrequency(find(~isnan(data_1D.AqFrequency), 1, 'first')))
    meanFrequency = 0;
    meanFreqStr = sprintf('mean offset to AQ.Frequency = %.2f Hz', meanOffset);
    ylabelStr = {'Frequency', 'in Hz'};
  else
    meanFrequency = meanNAN(data_1D.AqFrequency);
    meanFreqStr = sprintf('mean offset to %.6f MHz = %.2f Hz', data_1D.AqFrequency(1)/1e6, meanOffset);
    ylabelStr = {'f_{Offset}', 'in Hz'};
  end
  if numel(h_lines) ~= 1 || ~ishghandle(h_lines, 'line') || ...
      isempty(hMeanFreq) || ~ishghandle(hMeanFreq, 'textboxshape') || ...
      length(hLabelY) < 3 || ~ishghandle(hLabelY(3), 'text')
    % A former annotation is hidden somewhere in the figure descendents
    delete(findall(hFigure, 'Tag', 'plotData1DFreqAnnotation'));
    h_lines = plot(hAxes(3), ...
      data_1D.(config.timeFieldname), signalFrequency-meanFrequency, ...
      'Tag', 'plotData1D_LineFreq');
    grid(hAxes(3), 'on');
    isNewAxes = true;
    hLabelY(3) = ylabel(hAxes(3), ylabelStr);
    pos = get(hAxes(3), 'Position');
    hMeanFreq = annotation(hParent, 'textbox', ...
      'Position', [pos(1), pos(2)+0.8*pos(4), pos(3), 0.2*pos(4)], ...
      'String', meanFreqStr, 'HorizontalAlignment', 'center', 'EdgeColor', 'none', ...
      'VerticalAlignment', 'top', 'Margin', 2, 'Tag', 'plotData1DFreqAnnotation');
    set(hAxes(3), 'Tag', 'plotData1D_Axes3', 'DeleteFcn', {@DeleteFreqAxes, hMeanFreq});
  else
    set(h_lines, 'XData', data_1D.(config.timeFieldname), 'YData', signalFrequency-meanFrequency)
    set(hLabelY(3), 'String', ylabelStr)
    set(hMeanFreq, 'String', meanFreqStr);
  end
  if length(data_1D.(config.timeFieldname))/sum(isnan(data_1D.(config.timeFieldname)))==2 && ~config.omitnan
    set(h_lines, 'LineStyle', 'none', 'Marker', '.');
  else
    set(h_lines, 'LineStyle', '-', 'Marker', 'none');
  end
end

if numAxesChanged || isNewAxes
  validAxes = hAxes(plotAxes);
  set(validAxes(end), 'XTickLabelMode', 'auto');
  xlabel(validAxes(end), 'time in s');
  if numel(validAxes) > 1
    set(validAxes(1:end-1), 'XTickLabel', '');
    % re-arrange axes
    drawnow('expose'); % necessary to have updated axes positions below; FIXME: Can this be avoided?
    % Decrease vertical spacing between axes
    axesPos = get(validAxes(1), 'Position');
    top = axesPos(2) + axesPos(4);
    pos = get(validAxes(end), 'Position');
    bottom = pos(2);
    spacerPart = 0.25; % portion of axes height that is used as spacer between axes
    numAxes = numel(validAxes);
    axesHeight = (top-bottom) / (numAxes + (numAxes-1)*spacerPart);
    for iAxes = 1:numAxes
      axesPos(4) = axesHeight;
      axesPos(2) = (numAxes-iAxes) * (1+spacerPart)*axesHeight + bottom;
      set(validAxes(iAxes), 'Position', axesPos);
    end
  else
    subplot(1,1,1, validAxes);
  end
  if config.plotFreq
    pos = get(hAxes(3), 'Position');
    set(hMeanFreq, 'Position', [pos(1), pos(2)+0.8*pos(4), pos(3), 0.2*pos(4)]);
  end
end
% clear HW data_1D
% refresh(hParent)
% drawnow expose

end


function DeleteFreqAxes(hObject, eventData, hMeanFreq)
%% Delete annotation with frequency axes

if ishghandle(hMeanFreq)
  delete(hMeanFreq);
end

end


function outputTxt = plotData1D_DataCursorFcn(pointDataTip, eventObj, HW, iDevice)
% pointDataTip A PointDataTip object (undocumented Matlab and currently not used)
% eventObj     Object containing event data structure
% outputTxt    Data cursor text (string or cell array of strings)
pos = get(eventObj, 'Position');
target = get(eventObj, 'Target');
disp('t')
targetTag = get(target, 'Tag');
switch targetTag
  case {'plotData1D_Line1', 'plotData1D_Line2', 'plotData1D_Line3'}
    % data axes
    outputTxt = {['time: ' num2str(pos(1), 3) ' s'], ...
      ['amplitude: ' num2str(pos(2), 3) ' ' HW.RX(iDevice).AmplitudeUnit]};
  case 'plotData1D_LinePhase'
    % phase axes
    outputTxt = {['time: ' num2str(pos(1), 3) ' s'], ['phase: ' num2str(pos(2), 3) ' rad']};
  case 'plotData1D_LineFreq'
    % frequency (offset)
    ylab = get(get(get(target, 'Parent'), 'YLabel'), 'String');
    outputTxt = {['time: ' num2str(pos(1), 3) ' s'], [ylab{1} ': ' num2str(pos(2), 3) ' Hz']};
  otherwise
    outputTxt = [];
end

end


function plotData1D_DeleteFcn(hParent, eventData)
%% Executes on deletion of parent (panel or figure).
%
%       plotSeq_DeleteFcn(hParent, eventData)
%
% This function un-registers the DataCursorFcn when the graphics parent is
% deleted.

% FIXME: handle multiple plotData1D in one figure
hFigure = ancestor(hParent, 'figure');
DataCursorUpdateFcns = getappdata(hFigure, 'DataCursorUpdateFcns');
if isfield(DataCursorUpdateFcns, 'plotData1D')
  DataCursorUpdateFcns = rmfield(DataCursorUpdateFcns, 'plotData1D');
  setappdata(hFigure, 'DataCursorUpdateFcns', DataCursorUpdateFcns);
end

end
