function plotSequence(HW, Seq, AQ, TX, Grad)
%% Plot sequence with ui controls
%
%     plotSequence(HW, Seq, AQ, TX, Grad)
%
% This function opens a figure showing the pulse program (rf pulses, acquisition
% windows, gradient pulses and digital IO signals) of the sequence. For this,
% the function "plotSeq" is used. Additionally, some controls are displayed in
% the upper right corner of the parent figure or uipanel.
%
% INPUT:
%   See help for function plotSeq.
%   Additional setting in Seq.plotSequence or settings with differing defaults:
%     Gradients       Can be used to override the setting in Seq.plotSeq
%                     (Default: Seq.plotSeq).
%     showDuration    Show duration of sequence in side bar (default: true).
%     possibleWraps   Array with number of possible wraps. If this is
%                     non-scalar, a drop-down list is used to display these
%                     values instead of a text-box for wraps (default:
%                     Seq.plotSeq.wraps).
%     AQColors        1x2 cell array of color values for each AQ channel
%                     (default: the 4th & 5th color of "DefaultAxesColorOrder").
%     controls        Structure with the following fields (defaults in
%                     parentesis) that controls which control interfaces are
%                     displayed:
%       Gradients (1:3), stackGrads (true), stackTXRX (true), wraps (true),
%       tRepStart (true), tRepEnd (true)
%       figureModes     Display a button group on the side panel to control the
%                       figure mode ("pan", "zoom", "data cursor", or "none").
%                       (Default: HW.PlotSequence.controls.figureModes = false).
%
% ------------------------------------------------------------------------------
% (C) Copyright 2017-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% default input
Seq = set_EmptyField(Seq, 'plotSeq', 1:3);
Seq = set_EmptyField(Seq, 'plotSeqStart', 1);
Seq = set_EmptyField(Seq, 'plotSeqEnd', length(Seq.tRep));
Seq = set_EmptyField(Seq, 'externalLoops', 1);

% structure with configuration
Seq = set_EmptyField(Seq, 'plotSequence', HW.PlotSequence);
PlotSequence = Seq.plotSequence;
PlotSequence = set_EmptyField(PlotSequence, 'Gradients', Seq.plotSeq);
Seq.plotSeq = PlotSequence.Gradients;
PlotSequence = set_EmptyField(PlotSequence, 'hParent', HW.PlotSequence.hParent);
PlotSequence = set_EmptyField(PlotSequence, 'wraps', HW.PlotSequence.wraps);
PlotSequence = set_EmptyField(PlotSequence, 'possibleWraps', PlotSequence.wraps);
PlotSequence = set_EmptyField(PlotSequence, 'possibleWrapsNames', arrayfun(@num2str, PlotSequence.possibleWraps, 'UniformOutput', false));
PlotSequence = set_EmptyField(PlotSequence, 'stackGrads', HW.PlotSequence.stackGrads);
PlotSequence = set_EmptyField(PlotSequence, 'stackTXRX', HW.PlotSequence.stackTXRX);
PlotSequence = set_EmptyField(PlotSequence, 'raiseFigure', HW.PlotSequence.raiseFigure);

colorOrder = get(0, 'DefaultAxesColorOrder');
colorOrder = mat2cell(colorOrder, ones(size(colorOrder,1),1), 3);
PlotSequence = set_EmptyField(PlotSequence, 'AQColors', colorOrder(5:6));

PlotSequence = set_EmptyField(PlotSequence, 'controls', struct());
PlotSequence.controls = set_EmptyField(PlotSequence.controls, 'Gradients', 1:3);
PlotSequence.controls = set_EmptyField(PlotSequence.controls, 'stackGrads', true);
PlotSequence.controls = set_EmptyField(PlotSequence.controls, 'stackTXRX', true);
PlotSequence.controls = set_EmptyField(PlotSequence.controls, 'wraps', true);
PlotSequence.controls = set_EmptyField(PlotSequence.controls, 'tRepStart', true);
PlotSequence.controls = set_EmptyField(PlotSequence.controls, 'tRepEnd', true);
PlotSequence.controls = set_EmptyField(PlotSequence.controls, 'figureModes', HW.PlotSequence.controls.figureModes);
PlotSequence.controls = set_EmptyField(PlotSequence.controls, 'resetView', true);
PlotSequence = set_EmptyField(PlotSequence, 'showDuration', false);

hParent = PlotSequence.hParent;

%% open (and clear) parent figure if necessary
isFigure = false;
if ishghandle(hParent, 'figure') || (isa(hParent, 'double') && mod(hParent, 1) == 0)
  if PlotSequence.raiseFigure || ~ishghandle(hParent, 'figure')
    hFigure = figure(hParent);
  else
    hFigure = hParent;
  end
  isFigure = true;
  set(hFigure, 'Name', 'Pulse Program');
elseif ishghandle(hParent, 'uipanel')
  hFigure = ancestor(hParent, 'figure');
else
  error('PD:plotSequence', ...
    '"Seq.plotSequence.hParent" must be a valid handle to a figure or uipanel.');
end
oldPlotSequence = getappdata(hParent, 'plotSequence');
if isempty(oldPlotSequence), oldPlotSequence = PlotSequence; end

%% find uipanels
hPanels = getappdata(hParent, 'plotSequencePanels');

newControls = false;
if numel(hPanels) ~= 3 || ~all(ishghandle(hPanels, 'uipanel')) || ...
    oldPlotSequence.controls.figureModes ~= PlotSequence.controls.figureModes
  if isFigure
    clf(hFigure);
  else
    delete(get(hParent, 'Children'));
  end
  % use approximate positions (do not really matter yet)
  hPanels(1) = uipanel(hParent, 'Position', [0, 0, 0.8, 1], 'BorderType', 'none', 'Visible', 'off', 'Tag', 'plotSequenceMainPanel');
  hPanels(2) = uipanel(hParent, 'Position', [0.8, 0.2, 0.2, 0.8], 'DeleteFcn', @plotSequence_DeleteFcn, 'Visible', 'off');
  hPanels(3) = uipanel(hParent, 'Position', [0.8, 0.2, 0, 0.2], 'BorderType', 'none', 'Tag', 'ErrorMessages', 'Visible', 'off');
  set(hParent, 'ResizeFcn', @parentResizeFcn);
  newControls = true;
else
  % settings in GUI have priority over settings in sequence
  oldWrapIdx = find(oldPlotSequence.wraps == oldPlotSequence.possibleWraps, 1, 'first');
  if ~isempty(oldWrapIdx) && numel(PlotSequence.possibleWraps) >= oldWrapIdx
    newWraps = PlotSequence.possibleWraps(oldWrapIdx);
  else
    newWraps = PlotSequence.wraps;
  end
  newPossibleWraps = PlotSequence.possibleWraps;
  PlotSequence = oldPlotSequence;
  PlotSequence.wraps = newWraps;
  PlotSequence.possibleWraps = newPossibleWraps;
  % settings in dedicated structure have priority over settings in sequence
  Seq.plotSeq = PlotSequence.Gradients;
end
if isFigure
  set(hFigure, 'Name', 'PulseProgram');
end

%% main panel with Pulse Program
PlotSequence.hParent = hPanels(1);
Seq.plotSequence = PlotSequence;
setappdata(hParent, 'plotSequence', PlotSequence);
plotSeq(HW, Seq, AQ, TX, Grad);

%% side panel with controls
setappdata(hParent, 'plotSequencePanels', hPanels);

SeqData.HW = HW;
SeqData.Seq = Seq;
SeqData.AQ = AQ;
SeqData.TX = TX;
SeqData.Grad = Grad;
setappdata(hPanels(2), 'plotSequenceData', SeqData);

if newControls
  numLines = sum(PlotSequence.controls.Gradients>0 & PlotSequence.controls.Gradients<5) + ...
             PlotSequence.controls.stackGrads + ...
             PlotSequence.controls.stackTXRX + ...
             2 * (1 + PlotSequence.controls.tRepStart + PlotSequence.controls.tRepEnd) * PlotSequence.controls.wraps + ...
             PlotSequence.controls.resetView + ...
             2 * PlotSequence.showDuration;
  numLinesFigModesGroup = 4;
  if PlotSequence.controls.figureModes
    numLines = numLines + numLinesFigModesGroup;
  end
  setappdata(hPanels(2), 'plotSequenceNumLines', numLines);

  lineNum = 1;
  if any(PlotSequence.controls.Gradients == 1)
    hControls.plotGrad1 = uicontrol(hPanels(2), 'Units', 'normalized', ...
      'Position', [0, 1-lineNum/numLines, 1, 1/numLines], ...
      'Style', 'checkbox', 'String', HW.Grad.Name{1}, 'Value', any(PlotSequence.Gradients==1), ...
      'Callback', {@CallbackCheckboxGrad, 1}, 'Interruptible', 'off');
    lineNum = lineNum+1;
  end
  if any(PlotSequence.controls.Gradients == 2)
    hControls.plotGrad2 = uicontrol(hPanels(2), 'Units', 'normalized', ...
      'Position', [0, 1-lineNum/numLines, 1, 1/numLines], ...
      'Style', 'checkbox', 'String', HW.Grad.Name{2}, 'Value', any(PlotSequence.Gradients==2), ...
      'Callback', {@CallbackCheckboxGrad, 2}, 'Interruptible', 'off');
    lineNum = lineNum+1;
  end
  if any(PlotSequence.controls.Gradients == 3)
    hControls.plotGrad3 = uicontrol(hPanels(2), 'Units', 'normalized', ...
      'Position', [0, 1-lineNum/numLines, 1, 1/numLines], ...
      'Style', 'checkbox', 'String', HW.Grad.Name{3}, 'Value', any(PlotSequence.Gradients==3), ...
      'Callback', {@CallbackCheckboxGrad, 3}, 'Interruptible', 'off');
    lineNum = lineNum+1;
  end
  if any(PlotSequence.controls.Gradients == 4)
    hControls.plotGrad4 = uicontrol(hPanels(2), 'Units', 'normalized', ...
      'Position', [0, 1-lineNum/numLines, 1, 1/numLines], ...
      'Style', 'checkbox', 'String', HW.Grad.Name{4}, 'Value', any(PlotSequence.Gradients==4), ...
      'Callback', {@CallbackCheckboxGrad, 4}, 'Interruptible', 'off');
    lineNum = lineNum+1;
  end

  if PlotSequence.controls.stackGrads
    hControls.stackGrads = uicontrol(hPanels(2), 'Units', 'normalized', ...
      'Position', [0, 1-lineNum/numLines, 1, 1/numLines], ...
      'Style', 'checkbox', 'String', 'Stack Gradients', 'Value', PlotSequence.stackGrads, ...
      'Callback', {@CallbackCheckboxStack, 'stackGrads'}, 'Interruptible', 'off');
    lineNum = lineNum+1;
  end

  if PlotSequence.controls.stackTXRX
    hControls.stackTXRX = uicontrol(hPanels(2), 'Units', 'normalized', ...
      'Position', [0, 1-lineNum/numLines, 1, 1/numLines], ...
      'Style', 'checkbox', 'String', 'Stack TX/RX', 'Value', PlotSequence.stackTXRX, ...
      'Callback', {@CallbackCheckboxStack, 'stackTXRX'}, 'Interruptible', 'off');
    lineNum = lineNum+1;
  end

  if PlotSequence.controls.wraps
    if PlotSequence.controls.tRepStart
      uicontrol(hPanels(2), 'Units', 'normalized', ...
        'Position', [0, 1-lineNum/numLines, 1, 1/numLines], ...
        'Style', 'text', 'String', 'First tRep in plot:', ...
        'HorizontalAlignment', 'left');
      lineNum = lineNum+1;
      hControls.plotSeqStart = uicontrol(hPanels(2), 'Units', 'normalized', ...
        'Position', [0.1, 1-lineNum/numLines, 0.8, 1/numLines], ...
        'Style', 'edit', 'String', num2str(Seq.plotSeqStart, '%d'), ...
        'Callback', {@CallbackEdit, 'plotSeqStart'}, 'Interruptible', 'off');
      lineNum = lineNum+1;
    end

    if PlotSequence.controls.tRepEnd
      uicontrol(hPanels(2), 'Units', 'normalized', ...
        'Position', [0, 1-lineNum/numLines, 1, 1/numLines], ...
        'Style', 'text', 'String', 'Last tRep in plot:', ...
        'HorizontalAlignment', 'left');
      lineNum = lineNum+1;
      hControls.plotSeqEnd = uicontrol(hPanels(2), 'Units', 'normalized', ...
        'Position', [0.1, 1-lineNum/numLines, 0.8, 1/numLines], ...
        'Style', 'edit', 'String', num2str(Seq.plotSeqEnd, '%d'), ...
        'Callback', {@CallbackEdit, 'plotSeqEnd'}, 'Interruptible', 'off');
      lineNum = lineNum+1;
    end

    uicontrol(hPanels(2), 'Units', 'normalized', ...
      'Position', [0, 1-lineNum/numLines, 1, 1/numLines], ...
      'Style', 'text', 'String', 'Number of Wraps:', 'HorizontalAlignment', 'left');
    lineNum = lineNum+1;
    if numel(Seq.plotSequence.possibleWraps) > 1
      selectedWrap = find(Seq.plotSequence.wraps == Seq.plotSequence.possibleWraps, 1, 'first');
      if isempty(selectedWrap), selectedWrap = 1; end
      hControls.plotSeqWraps = uicontrol(hPanels(2), 'Units', 'normalized', ...
        'Position', [0.1, 1-(lineNum-1)/numLines-0.6/numLines, 0.8, 0.8/numLines], ...
        'Style', 'popupmenu', 'Value', selectedWrap, ...
        'String', Seq.plotSequence.possibleWrapsNames, ...
        'Callback', {@CallbackPopupmenu, 'wraps', 'possibleWraps'}, 'Interruptible', 'off');
    else
      hControls.plotSeqWraps = uicontrol(hPanels(2), 'Units', 'normalized', ...
        'Position', [0.1, 1-lineNum/numLines, 0.8, 1/numLines], ...
        'Style', 'edit', 'String', num2str(Seq.plotSequence.wraps, '%d'), ...
        'Callback', {@CallbackEdit, 'wraps'}, 'Interruptible', 'off');
    end
    lineNum = lineNum+1;
  end

  % buttongroup
  if PlotSequence.controls.figureModes
    lineNum = lineNum-1+numLinesFigModesGroup;
    hControls.modeControls = uibuttongroup(hPanels(2), 'Units', 'normalized', ...
      'Position', [0, 1-lineNum/numLines, 1, numLinesFigModesGroup/numLines], ...
      'SelectionChangedFcn', @CallbackFigureModeChanged);
    lineNumGroup = 1;
    hControls.modes.None = uicontrol(hControls.modeControls, 'Units', 'normalized', ...
      'Position', [0, 1-lineNumGroup/numLinesFigModesGroup, 1, 1/numLinesFigModesGroup], ...
      'Style', 'radiobutton', 'String', 'none', 'Tag', 'none');
    lineNumGroup = lineNumGroup+1;
    hControls.modes.Zoom = uicontrol(hControls.modeControls, 'Units', 'normalized', ...
      'Position', [0, 1-lineNumGroup/numLinesFigModesGroup, 1, 1/numLinesFigModesGroup], ...
      'Style', 'radiobutton', 'String', 'zoom', 'Tag', 'zoom');
    lineNumGroup = lineNumGroup+1;
    hControls.modes.Pan = uicontrol(hControls.modeControls, 'Units', 'normalized', ...
      'Position', [0, 1-lineNumGroup/numLinesFigModesGroup, 1, 1/numLinesFigModesGroup], ...
      'Style', 'radiobutton', 'String', 'pan', 'Tag', 'pan');
    lineNumGroup = lineNumGroup+1;
    hControls.modes.DataCursor = uicontrol(hControls.modeControls, 'Units', 'normalized', ...
      'Position', [0, 1-lineNumGroup/numLinesFigModesGroup, 1, 1/numLinesFigModesGroup], ...
      'Style', 'radiobutton', 'String', 'data cursor', 'Tag', 'datacursor');
    lineNum = lineNum + 1;
  end

  if PlotSequence.controls.figureModes
    % check current figure mode
    hzm = zoom(hFigure);
    hpm = pan(hFigure);
    hdcm = datacursormode(hFigure);

    if strcmp(hzm.Enable, 'on')
      set(hControls.modes.Zoom, 'Value', 1);
    elseif strcmp(hpm.Enable, 'on')
      set(hControls.modes.Pan, 'Value', 1);
    elseif strcmp(hdcm.Enable, 'on')
      set(hControls.modes.DataCursor, 'Value', 1);
    else
      set(hControls.modes.None, 'Value', 1);
    end
  end
  % figure mode listener
  hManager = uigetmodemanager(hFigure);
  fml = addlistener(hManager, 'CurrentMode', 'PostSet', ...
    @(hObj,evt) figureModeListener(hObj, evt, hControls.modes));
  setappdata(hPanels(2), 'plotSequenceFigureModeListener', fml);

  if PlotSequence.controls.resetView
    hControls.resetView = uicontrol(hPanels(2), 'Units', 'normalized', ...
      'Position', [0, 1-lineNum/numLines, 1, 1/numLines], ...
      'Style', 'pushbutton', ...
      'String', 'Reset View', ...
      'Callback', @CallbackResetView);
    % lineNum = lineNum+1;
  else
    hControls.resetView = [];
  end

  % display sequence duration
  if PlotSequence.showDuration
    uicontrol(hPanels(2), 'Units', 'normalized', ...
      'Position', [0, 1-lineNum/numLines, 1, 1/numLines], ...
      'Style', 'text', 'String', 'Duration of sequence:', ...
      'HorizontalAlignment', 'left');
    lineNum = lineNum + 1;
    hControls.duration = uicontrol(hPanels(2), 'Units', 'normalized', ...
      'Position', [0, 1-lineNum/numLines, 1, 1/numLines], ...
      'Style', 'text', 'String', sprintf('%.3f s', (sum(Seq.tRep(:))+Seq.TimeToNextSequence)*Seq.externalLoops), ...
      'HorizontalAlignment', 'left');
  end
  setappdata(hPanels(2), 'plotSequenceControls', hControls);

else
  hControls = getappdata(hPanels(2), 'plotSequenceControls');
  SetIfValid(hControls, 'plotGrad1', 'Value', any(Seq.plotSeq==1));
  SetIfValid(hControls, 'plotGrad2', 'Value', any(Seq.plotSeq==2));
  SetIfValid(hControls, 'plotGrad3', 'Value', any(Seq.plotSeq==3));
  SetIfValid(hControls, 'plotGrad4', 'Value', any(Seq.plotSeq==4));
  SetIfValid(hControls, 'stackGrads', 'Value', PlotSequence.stackGrads);
  SetIfValid(hControls, 'stackTXRX', 'Value', PlotSequence.stackTXRX);
  SetIfValid(hControls, 'plotSeqStart', 'String', num2str(Seq.plotSeqStart, '%d'));
  SetIfValid(hControls, 'plotSeqEnd', 'String', num2str(Seq.plotSeqEnd, '%d'));
  if PlotSequence.controls.wraps && numel(Seq.plotSequence.possibleWraps) > 1
    set(hControls.plotSeqWraps, 'String', arrayfun(@num2str, Seq.plotSequence.possibleWraps, 'UniformOutput', false));
    oldIdx = get(hControls.plotSeqWraps, 'Value');
    if oldIdx > numel(Seq.plotSequence.possibleWraps)
      newIdx = find(Seq.plotSequence.possibleWraps == Seq.plotSequence.wraps, 1, 'first');
      if isempty(newIdx), newIdx = 1; end
    else
      newIdx = oldIdx;
    end
    SetIfValid(hControls, 'plotSeqWraps', ...
      'String', Seq.plotSequence.possibleWrapsNames, ...
      'Value', newIdx);
  else
    SetIfValid(hControls, 'plotSeqWraps', 'String', num2str(Seq.plotSequence.wraps, '%d'));
  end
  SetIfValid(hControls, 'duration', 'String', sprintf('%.3f s', (sum(Seq.tRep(:))+Seq.TimeToNextSequence)*Seq.externalLoops));
end


%% finalize parent
parentResizeFcn(hParent);
set(hPanels, 'Visible', 'on');

end

function SetIfValid(inStr, field, varargin)
%% set inStr.(field) if it is valid
if isfield(inStr, field) && ishghandle(inStr.(field))
  set(inStr.(field), varargin{:});
end

end


function CallbackCheckboxGrad(hObject, eventData, numGrad)
%% Update gradients

SeqData = getappdata(get(hObject, 'Parent'), 'plotSequenceData');
plotSequence = getappdata(get(get(hObject, 'Parent'), 'Parent'), 'plotSequence');

if get(hObject, 'Value')
  SeqData.Seq.plotSeq = unique([SeqData.Seq.plotSeq,numGrad]);
  SeqData.Seq.plotSeq(SeqData.Seq.plotSeq==0) = [];
  plotSequence.Gradients = unique([plotSequence.Gradients,numGrad]);
  plotSequence.Gradients(plotSequence.Gradients==0) = [];
else
  SeqData.Seq.plotSeq(SeqData.Seq.plotSeq == numGrad) = [];
  if isempty(SeqData.Seq.plotSeq)
    SeqData.Seq.plotSeq = 0;
  end
  plotSequence.Gradients(plotSequence.Gradients == numGrad) = [];
  if isempty(plotSequence.Gradients)
    plotSequence.Gradients = 0;
  end
end
SeqData.Seq.plotSequence = plotSequence;
setappdata(get(hObject, 'Parent'), 'plotSequenceData', SeqData);
setappdata(get(get(hObject, 'Parent'), 'Parent'), 'plotSequence', plotSequence);

plotSeq(SeqData.HW, SeqData.Seq, SeqData.AQ, SeqData.TX, SeqData.Grad);

end


function CallbackCheckboxStack(hObject, eventData, type)
%% Update stack callbacks
SeqData = getappdata(get(hObject, 'Parent'), 'plotSequenceData');
plotSequence = getappdata(get(get(hObject, 'Parent'), 'Parent'), 'plotSequence');

plotSequence.(type) = get(hObject, 'Value');
SeqData.Seq.plotSequence = plotSequence;
setappdata(get(hObject, 'Parent'), 'plotSequenceData', SeqData);
setappdata(get(get(hObject, 'Parent'), 'Parent'), 'plotSequence', plotSequence);

plotSeq(SeqData.HW, SeqData.Seq, SeqData.AQ, SeqData.TX, SeqData.Grad);

end


function CallbackEdit(hObject, eventData, type)
%% Update edit boxes

SeqData = getappdata(get(hObject, 'Parent'), 'plotSequenceData');

newVal = str2double(get(hObject, 'String'));
if isnan(newVal)
  switch type
    case 'plotSeqStart'
      newVal = SeqData.Seq.plotSeqStart;
    case 'plotSeqEnd'
      newVal = SeqData.Seq.plotSeqEnd;
    case 'wraps'
      newVal = SeqData.Seq.plotSequence.wraps;
  end
  set(hObject, 'String', newVal);
  return;
end
switch type
  case 'plotSeqStart'
    SeqData.Seq.plotSeqStart = newVal;
  case 'plotSeqEnd'
    SeqData.Seq.plotSeqEnd = newVal;
  case 'wraps'
    SeqData.Seq.plotSequence.wraps = newVal;
end
setappdata(get(hObject, 'Parent'), 'plotSequenceData', SeqData);
setappdata(get(get(hObject, 'Parent'), 'Parent'), 'plotSequence', SeqData.Seq.plotSequence);

plotSeq(SeqData.HW, SeqData.Seq, SeqData.AQ, SeqData.TX, SeqData.Grad);

end


function CallbackPopupmenu(hObject, eventData, type, valuesField)
%% Update edit boxes

SeqData = getappdata(get(hObject, 'Parent'), 'plotSequenceData');

newIdx = get(hObject, 'Value');

if ~isfield(SeqData.Seq.plotSequence, valuesField)
  errordlg('Error in plotSequence::CallbackPopupmenu.');
  error('No such field %s', valuesField);
end

possibleValues = SeqData.Seq.plotSequence.(valuesField);
newVal = possibleValues(newIdx);

switch type
  case 'wraps'
    SeqData.Seq.plotSequence.wraps = newVal;
end
setappdata(get(hObject, 'Parent'), 'plotSequenceData', SeqData);
setappdata(get(get(hObject, 'Parent'), 'Parent'), 'plotSequence', SeqData.Seq.plotSequence);

plotSeq(SeqData.HW, SeqData.Seq, SeqData.AQ, SeqData.TX, SeqData.Grad);

end


function CallbackFigureModeChanged(hObject, eventData)
%% Change figure mode

hFigure = ancestor(hObject, 'figure');

switch get(eventData.NewValue, 'Tag')
  case 'zoom'
    zoom(hFigure, 'on');

  case 'pan'
    pan(hFigure, 'on');

  case 'datacursor'
    hdcm = datacursormode(hFigure);
    set(hdcm, 'Enable', 'on');
    % Re-install data cursor update function. It might be lost if the figure
    % ancestor has changed.
    % FIXME: Can we track such a change? Can we un-register from the old figure?
    HW = PD.HW.GetInstance();  % FIXME: This only works if HW is an object.
    DataCursorUpdateFcns = getappdata(hFigure, 'DataCursorUpdateFcns');
    % FIXME: How to get the correct iDevice?
    iDevice = 1;
    DataCursorUpdateFcns.plotSeq = @(pointDataTip, eventData) plotSeq_DataCursorFcn(pointDataTip, eventData, HW, iDevice);
    setappdata(hFigure, 'DataCursorUpdateFcns', DataCursorUpdateFcns);
    set(hdcm, 'UpdateFcn', @DataCursorUpdateFcnHandler);

  otherwise
    switch eventData.OldValue.String
      case 'zoom'
        zoom(hFigure, 'off');
      case 'pan'
        pan(hFigure, 'off');
      case 'data cursor'
        datacursormode(hFigure, 'off');
    end

end

end


function CallbackResetView(hObject, eventData)
%% Reset the view of the axes in the main panel

% get handle to the main panel
hParent = get(get(hObject, 'Parent'), 'Parent');
hPanels = getappdata(hParent, 'plotSequencePanels');

% find axes in main panel
hAxes = findobj(get(hPanels(1), 'Children'), 'flat', 'Type', 'axes');

% FIXME: This is an un-documented function.
resetplotview(hAxes, 'ApplyStoredView');

end


function parentResizeFcn(hObject, eventData)
%% Resize function

parentUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
parentPosition = get(hObject, 'Position');
set(hObject, 'Units', parentUnits);

hPanels = getappdata(hObject, 'plotSequencePanels');

if numel(hPanels) < 3
  return;
end

controlWidth = 120;
controlHeight = getappdata(hPanels(2), 'plotSequenceNumLines')*18; % line height: 18px
minMainWidth = 150;
betterControlWidth = max(0, min(controlWidth, parentPosition(3)-minMainWidth));
set(hPanels(1), 'Units', 'pixels', ...
  'Position', [0, 0, parentPosition(3) - betterControlWidth, parentPosition(4)]);
set(hPanels(2), 'Units', 'pixels', ...
  'Position', [parentPosition(3)-betterControlWidth, parentPosition(4)-controlHeight, ...
               betterControlWidth, controlHeight]);

set(hPanels(3), 'Units', 'pixels');
oldPos = get(hPanels(3), 'Position');
set(hPanels(3), 'Units', 'pixels', ...
  'Position', [parentPosition(3)-betterControlWidth, 0, ...
               betterControlWidth, oldPos(4)]);
end


function plotSequence_DeleteFcn(hObject, eventData)
%% Remove ResizeFcn from parent

hParent = get(hObject, 'Parent');
set(hParent, 'ResizeFcn', '');
hFigure = ancestor(hObject, 'figure');
zoom(hFigure, 'off');
pan(hFigure, 'off');
datacursormode(hFigure, 'off');

fml = getappdata(hObject, 'plotSequenceFigureModeListener');
if ~isempty(fml)
  delete(fml);
end

end


function figureModeListener(hObject, eventData, modeControls)
%% Listener to figure mode changes

persistent oldMode

hManager = eventData.AffectedObject;

if ~isempty(hManager.CurrentMode) && ...
    ~isequal(ancestor(modeControls.Zoom, 'figure'), hManager.CurrentMode.FigureHandle)
  % Re-install figure mode listener. It might be attached to the wrong figure if
  % the figure ancestor has changed.
  % FIXME: Can we track such a change?

  % remove old figure mode listener
  hControlsPanel = get(get(modeControls.Zoom, 'Parent'), 'Parent');
  fml = getappdata(hControlsPanel, 'plotSequenceFigureModeListener');
  delete(fml);
  % install new figure mode listener
  fml = addlistener(hManager, 'CurrentMode', 'PostSet', ...
    @(hObj,evt) figureModeListener(hObj, evt, modeControls));
  setappdata(hControlsPanel, 'plotSequenceFigureModeListener', fml);
  return;
end

if isempty(hManager.CurrentMode)
  set(modeControls.None, 'Value', 1);
  return;
end
switch hManager.CurrentMode.Name
  case 'Exploration.Zoom'
    set(modeControls.Zoom, 'Value', 1);
  case 'Exploration.Pan'
    set(modeControls.Pan, 'Value', 1);
  case 'Exploration.Datacursor'
    set(modeControls.DataCursor, 'Value', 1);
  otherwise
    set(modeControls.None, 'Value', 0);
    set(modeControls.Zoom, 'Value', 0);
    set(modeControls.Pan, 'Value', 0);
    set(modeControls.DataCursor, 'Value', 0);
end

if ~strcmp(oldMode, hManager.CurrentMode.Name)
  % FIXME: There is a bug in Matlab where PointDataTip objects become invalid when
  % the figure mode is changed to zoom which ultimately leads to a hard Matlab
  % crash.
  % Work around this by manually deleting all Data Cursors.
  if ~strcmp(hManager.CurrentMode.Name, 'Exploration.Datacursor')
    dataPoints = findall(ancestor(modeControls.Zoom, 'figure'), 'Type', 'hggroup');
    dataPoints(~arrayfun(@(x) isa(x, 'matlab.graphics.shape.internal.PointDataTip'), dataPoints)) = [];
    delete(dataPoints);
  end
end
oldMode = hManager.CurrentMode.Name;

end
