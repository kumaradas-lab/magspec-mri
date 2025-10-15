function plot_data_store(hf, plotData, hgProps)
%% Plot lines in subplots extending previous data
%
%       plot_data_store(hf, plotData, hgProps)
%
% INPUT:
%   hf          Handle to the figure or uipanel with the axes.
%
%   plotData    cell array with the data to plot. Every cell element corresponds
%               to one subplot. If the data is a vector, several lines are drawn
%               in the corresponding axes.
%
%   hgProps     structure with graphics properties of the following form
%               (all fields are optional):
%       hgProps.figure.props
%                   Properties that are directly passed to the figure. If hf is
%                   an uipanel, these settings are ignored.
%       hgProps.allAxes.props
%                   Properties that are directly passed to all axes.
%       hgProps.axes
%                   cell array of structs with the following fields that
%                   correspond to the subplots:
%       hgProps.axes{i}.props
%                       Properties that are directly passed to the
%                       respective axes.
%       hgProps.axes{i}.title
%                       String with the title of the respective subplot.
%       hgProps.axes{i}.xlabel
%                       String with the xlabel of the respective subplot.
%       hgProps.axes{i}.ylabel
%                       String with the ylabel of the respective subplot.
%       hgProps.axes{i}.legend.labels
%                       Cell array of strings with the legend labels for the
%                       respective subplot.
%       hgProps.axes{i}.legend.props
%                       Structure with the properties that are directly passed
%                       to the respective legend as property-value pairs.
%
% OUTPUT:
%   none
%
% To clear the stored data and to re-start with empty axes, call:
%       plot_data_store('clear')
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%%
persistent handles plotDataStore

%% early clear and return
if nargin == 1
    clear handles plotDataStore
    return
end

%% initialization
numAxes = numel(plotData);
[nRows, nCols] = size(plotData);

% add new plotData to plotDataStore
if isempty(plotDataStore) || ~iscell(plotDataStore) || ...
        numel(plotDataStore) ~= numAxes
    plotDataStore = cell(nRows, nCols);
end
for iAxes = 1:numAxes
    plotDataStore{iAxes}(:,end+1) = plotData{iAxes};
end

% adapt xlim if more iterations than set
extentedXLim = false;
if ~isemptyfield(hgProps, 'allAxes') && ...
        ~isemptyfield(hgProps.allAxes, 'props') && ...
        ~isempty(hgProps.allAxes.props)
    xlimFName = 'xlim';
    if isstruct(hgProps.allAxes.props)
        % find fieldname ignoring case
        fnames = fieldnames(hgProps.allAxes.props);
        xlimIdx = find(strcmpi(fnames, 'xlim'), 1, 'first');
        if ~isempty(xlimIdx)
            xlimFName = fnames{xlimIdx};
        end
    end
    if isfield(hgProps.allAxes.props, xlimFName) && ...
            hgProps.allAxes.props.(xlimFName)(2) < size(plotDataStore{iAxes}, 2)
        xlims = hgProps.allAxes.props.(xlimFName);
        xlims(2) = size(plotDataStore{1}, 2);
        extentedXLim = true;
    end
end

%% create figure with axes
if ~isfield(handles, 'axesPlot') || ...
    numel(handles.axesPlot) ~= numAxes || ...
    any(~ishghandle(handles.axesPlot, 'axes'))
  % find and clear or create figure
  if ishghandle(hf, 'figure')
    hfDataStore = clf(hf);
  elseif ~ishghandle(hf)
    % (re-)open closed figure
    hfDataStore = figure(hf);
  else
    % assume that "hf" is a valid parent
    hfDataStore = hf;
  end
  if ishghandle(hfDataStore, 'figure') && ...
      ~isemptyfield(hgProps, 'figure') && ...
      ~isemptyfield(hgProps.figure, 'props')
    set(hfDataStore, hgProps.figure.props);
  end

  % create axes
  for iAxes = 1:numAxes
    handles.axesPlot(iAxes) = subplot(nRows, nCols, iAxes, ...
      'Parent', hfDataStore);
    grid(handles.axesPlot(iAxes), 'on');
    if ~isemptyfield(hgProps, 'axes') && ...
        ~isemptyfield(hgProps.axes{iAxes}, 'title')
      title(handles.axesPlot(iAxes), hgProps.axes{iAxes}.title, ...
        'HandleVisibility', 'off');
    end
    if ~isemptyfield(hgProps, 'axes') && iscell(hgProps.axes) && ...
        numel(hgProps.axes)>=iAxes && ...
        ~isemptyfield(hgProps.axes{iAxes}, 'xlabel')
      xlabel(handles.axesPlot(iAxes), hgProps.axes{iAxes}.xlabel, ...
        'HandleVisibility', 'off');
    end
    if ~isemptyfield(hgProps, 'axes') && iscell(hgProps.axes) && ...
        numel(hgProps.axes)>=iAxes && ...
        ~isemptyfield(hgProps.axes{iAxes}, 'ylabel')
      ylabel(handles.axesPlot(iAxes), hgProps.axes{iAxes}.ylabel, ...
        'HandleVisibility', 'off');
    end
    if ~isemptyfield(hgProps, 'axes') && ...
        ~isemptyfield(hgProps.axes{iAxes}, 'props')
      set(handles.axesPlot(iAxes), hgProps.axes{iAxes}.props);
    end
  end
  % settings for all axes
  if ~isemptyfield(hgProps, 'allAxes') && ...
      ~isemptyfield(hgProps.allAxes, 'props')
    set(handles.axesPlot, hgProps.allAxes.props);
  end
  set(handles.axesPlot, 'NextPlot', 'replacechildren');
  linkaxes(handles.axesPlot, 'x');
end

%% plot in axes
for iAxes = 1:numAxes
  % find handles to lines in axes
  if ~isfield(handles, 'lines') || numel(handles.lines) < iAxes || ...
      ~all(ishghandle(handles.lines{iAxes}, 'line'))
    handles.lines{iAxes} = get(handles.axesPlot(iAxes), 'Children');
    if ~isempty(handles.lines{iAxes})
      handles.lines{iAxes}(~ishghandle(handles.lines{iAxes}, 'line')) = [];
    end
  end
  % plot in axes or update lines
  if isvector(plotDataStore{iAxes})
    if ~isscalar(handles.lines{iAxes})
      handles.lines{iAxes} = plot(handles.axesPlot(iAxes), plotDataStore{iAxes}');
    else
      set(handles.lines{iAxes}, 'YData', plotDataStore{iAxes})
    end
  else
    if numel(handles.lines{iAxes}) ~= size(plotDataStore{iAxes}, 1)
      handles.lines{iAxes} = plot(handles.axesPlot(iAxes), plotDataStore{iAxes}');
    else
      set(handles.lines{iAxes}, {'YData'}, mat2cell(plotDataStore{iAxes}, ...
        ones(size(handles.lines{iAxes})), size(plotDataStore{iAxes}, 2)))
    end
  end

  if extentedXLim
    set(handles.axesPlot(iAxes), 'XLim', xlims);
  end

  % legend
  if ~isemptyfield(hgProps.axes{iAxes}, 'legend') && ...
      (isemptyfield(handles, 'legends') || numel(handles.legends) < iAxes || ...
      ~ishghandle(handles.legends{iAxes}) || ...
      numel(get(handles.legends{iAxes}, 'String')) < numel(hgProps.axes{iAxes}.legend.labels))
    legendProps = fieldnames(hgProps.axes{iAxes}.legend.props).';
    legendProps(2,:) = cellfun(@(x) hgProps.axes{iAxes}.legend.props.(x), legendProps, 'UniformOutput', false);
    ws = warning('off', 'MATLAB:legend:IgnoringExtraEntries');
    handles.legends{iAxes} = legend(handles.axesPlot(iAxes), hgProps.axes{iAxes}.legend.labels, ...
      legendProps{:});
    if ~verLessThan('Matlab', '9.2')
      set(handles.legends{iAxes}, 'AutoUpdate', 'off'); % Do not automatically add new entries to legend
    end
    warning(ws);
  end
end

end
