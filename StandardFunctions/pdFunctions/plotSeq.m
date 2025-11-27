function plotSeq(HW, Seq, AQ, TX, Grad)
%% Plot sequence
%
%     plotSeq(HW, Seq, AQ, TX, Grad)
%
% This function opens a figure showing the pulse program (rf pulses, acquisition
% windows, gradient pulses and digital IO signals) of the sequence. The "tRep"s
% are shown consecutively starting at Seq.plotSeqStart and ending with
% Seq.plotSeqEnd unless "Seq.plotSequence.wraps" is other than 1 (see below).
%
%
% INPUT:
%
%   HW    HW structure or object
%
%   Seq   Structure with the following fields. If empty, default values are
%         used. These default values can be overwritten in HW.PlotSequence:
%     plotSeq         Array with the numbers of Gradients that are included in
%                     the pulse sequence plot (default: 1:3).
%                     Examples:
%                         Seq.plotSeq = 1:3;   % include x-, y-, and z-gradients
%                         Seq.plotSeq = [1,3]; % include x- and z-gradients
%                         Seq.plotSeq = 0;     % don't include any gradients
%                         Seq.plotSeq = [];    % Completely switch off automatic
%                                              % pulse sequence plot
%     plotSeqStart    Number of first tRep to include in the plot (default: 1).
%     plotSeqEnd      Number of last tRep to include in the plot (default:
%                     length(Seq.tRep) or such that the last wrap is complete).
%     plotSequence    structure with the following fields:
%         hParent     Handle to a parent (figure or uipanel) for the sequence
%                     plot. If empty, figure 21 is used.
%         figureName  If "hParent" is a figure, its title is set to this string
%                     (default: 'Pulse Program').
%         raiseFigure If false and the figure already exists, the figure is not
%                     raised above other windows and doesn't steal the focus
%                     when plotted (default: true).
%         wraps       Number of wraps in the plot (default: 1 = no wraps).
%         xLim        Limits for the x-axis (default: show complete sequence)
%         Gradients   Vector with the numbers of the gradient to include in the
%                     plot where "1" corresponds to the x-gradient, "2" to the y
%                     and "3" to the z-gradient. If empty or 0, the gradient
%                     pulses are not included in the plot. (Default: Seq.plotSeq
%                     which defaults to 1:3).
%         includeShim Add shim values to gradient plot (default: false).
%         stackGrads  If false, each gradient used its own axes. Otherwise, all
%                     gradients are plotted in the same axes (default: false).
%         stackTXRX   Use the same axes for the rf pulses and acquisition
%                     windows (default: false).
%         showTXPhase If false, only the amplitude of the TX pulse is plotted.
%                     If true, the phase of the pulse is indicated by its real
%                     and imaginary parts (default: true).
%         plotTX      Vector with indices of the TX structure that should be
%                     displayed (default: 1:numel(TX)).
%         plotDigitalIO
%                     Include digital IO in pulse program plot if there is
%                     activity at any channel
%                     (default: HW.PlotSequence.plotDigitalIO).
%         GradColors  1x4 cell array of color values for each gradient channel.
%                     If one cell contains several color values, these are
%                     linearly interpolated for the "stacked" tReps when using
%                     "wraps". When "wraps" is 1, the average color value is
%                     used. Default: x-gradients are reddish, y-gradients are
%                     greenish, z-gradients are blueish, and the fourth channel
%                     is of violetish hue.
%         TXColors    2x3 cell array of color values for each TX channel and the
%                     real, imaginary and absolute parts, respectively (default:
%                     the first six colors of "DefaultAxesColorOrder").
%         AQColors    1x2 cell array of color values for each AQ channel
%                     (default: {'b','r'}).
%         DigitalIOColor
%                     Color value for the Digital IO channels
%                     (default: the 7th color of "DefaultAxesColorOrder").
%         FontSize    Font size for axes labels and ticks. If empty, the default
%                     font sizes are used (default: []).
%         tOffset     Move all pulse program components in the graphics by an
%                     offset in seconds (default: 0).
%
%   AQ    structure with the settings for the acquisition windows
%
%   TX    structure with the settings for the rf pulses
%
%   Grad  structure with the settings for the gradients
%
%
% OUTPUT:
%   none
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
Seq = set_EmptyField(Seq, 'plotSeq', 1:3);
Seq = set_EmptyField(Seq, 'plotSeqStart', 1);
if Seq.plotSeqStart < 1, Seq.plotSeqStart = 1; end
Seq = set_EmptyField(Seq, 'plotSeqEnd', length(Seq.tRep));
if Seq.plotSeqEnd > length(Seq.tRep), Seq.plotSeqEnd = length(Seq.tRep); end

% structure with configuration
Seq = set_EmptyField(Seq, 'plotSequence', HW.PlotSequence);
plotSequence = Seq.plotSequence;
plotSequence = set_EmptyField(plotSequence, 'hParent', HW.PlotSequence.hParent);
plotSequence = set_EmptyField(plotSequence, 'wraps', HW.PlotSequence.wraps);
if plotSequence.wraps > Seq.plotSeqEnd-Seq.plotSeqStart+1
  plotSequence.wraps = Seq.plotSeqEnd-Seq.plotSeqStart+1;
end
if plotSequence.wraps < 1, plotSequence.wraps = 1; end
plotSequence = set_EmptyField(plotSequence, 'raiseFigure', HW.PlotSequence.raiseFigure);
plotSequence = set_EmptyField(plotSequence, 'stackGrads', HW.PlotSequence.stackGrads);
plotSequence = set_EmptyField(plotSequence, 'stackTXRX', HW.PlotSequence.stackTXRX);
plotSequence = set_EmptyField(plotSequence, 'Gradients', Seq.plotSeq);
plotSequence = set_EmptyField(plotSequence, 'figureName', 'Pulse Program');
plotSequence = set_EmptyField(plotSequence, 'includeShim', false);
plotSequence = set_EmptyField(plotSequence, 'showTXPhase', true);
if isemptyfield(plotSequence, 'plotTX')
  plotSequence.plotTX = 1:numel(TX);
end
plotSequence = set_EmptyField(plotSequence, 'plotDigitalIO', HW.PlotSequence.plotDigitalIO);
if ~isfield(plotSequence, 'FontSize'), plotSequence.FontSize = []; end
if isemptyfield(plotSequence, 'tOffset'), plotSequence.tOffset = 0; end

Seq.plotSeqEnd = floor((Seq.plotSeqEnd-Seq.plotSeqStart+1)/plotSequence.wraps)*plotSequence.wraps+Seq.plotSeqStart-1;
hParent = plotSequence.hParent;

% colors
plotSequence = set_EmptyField(plotSequence, 'GradColors', HW.PlotSequence.GradColors);
GradColors = plotSequence.GradColors;
for iColors = 1:length(GradColors)
  if iscellstr(GradColors{iColors})
    GradColors{iColors} = str2rgb(GradColors{iColors});
  end
  if iscell(GradColors{iColors})
    GradColors{iColors} = cell2mat(reshape(GradColors{iColors}, [], 1));
  end
end
plotSequence = set_EmptyField(plotSequence, 'TXColors', HW.PlotSequence.TXColors);
if isempty(plotSequence.TXColors) || isempty(plotSequence.DigitalIOColor)
%   colorOrder = get(0, 'DefaultAxesColorOrder');
  colorOrder =  [    0, 0.4470, 0.7410; ...
                0.8500, 0.3250, 0.0980; ...
                0.9290, 0.6940, 0.1250; ...
                0.4940, 0.1840, 0.5560; ...
                0.4660, 0.6740, 0.1880; ...
                0.3010, 0.7450, 0.9330; ...
                0.6350, 0.0780, 0.1840];
  colorOrder = mat2cell(colorOrder, ones(size(colorOrder,1),1), 3);
end
if isempty(plotSequence.TXColors)
  % use default colors for TX (abs, real, imag)
  plotSequence.TXColors = reshape(colorOrder(1:6), 3, 2).';
else
  plotSequence.TXColors = str2rgb(plotSequence.TXColors);
end
plotSequence = set_EmptyField(plotSequence, 'AQColors', HW.PlotSequence.AQColors);
plotSequence.AQColors = str2rgb(plotSequence.AQColors);
plotSequence = set_EmptyField(plotSequence, 'DigitalIOColor', HW.PlotSequence.DigitalIOColor);
if isempty(plotSequence.DigitalIOColor)
  % use default colors for DigitalIO
  plotSequence.DigitalIOColor = colorOrder{7};
else
  plotSequence.DigitalIOColor = str2rgb(plotSequence.DigitalIOColor);
end

if isemptyfield(plotSequence, 'iDevice'), plotSequence.iDevice = 1:numel(HW.MMRT); end

if isemptyfield(Seq, 'tOffset')
  Seq.tOffset = zeros(size(Seq.tRep));
end


%% handle multiple data channels
if ~isscalar(plotSequence.iDevice)
  plotSequenceIn = plotSequence;
  GradIn = Grad;
  for iDevice = plotSequence.iDevice(:).'
    if ishghandle(plotSequence.hParent, 'figure') || isnumeric(plotSequence.hParent)
      plotSequence.hParent = double(plotSequence.hParent) + 10*(iDevice-1);
    else
      plotSequence.hParent = HW.PlotSequence.hParent + 10*(iDevice-1);
    end
    plotSequence.iDevice = iDevice;
    plotSequence.figureName = sprintf('%s (device #%d)', plotSequence.figureName, iDevice);
    Seq.plotSequence = plotSequence;
    nGradDevices = [HW.Grad(:).n];
    nGradDevicesCum = [0, cumsum(nGradDevices)];
    Grad = GradIn((nGradDevicesCum(iDevice)+1):min(nGradDevicesCum(iDevice+1), numel(GradIn)));
    plotSeq(HW, Seq, AQ, TX, Grad);
    plotSequence = plotSequenceIn;
  end
  return;
end


%% open parent figure if necessary
if ishghandle(hParent, 'figure') || (isa(hParent, 'double') && mod(hParent, 1) == 0)
  if plotSequence.raiseFigure || ~ishghandle(hParent, 'figure')
    hFigure = figure(hParent);
  else
    hFigure = hParent;
  end
  if Seq.average>1
    set(hFigure, 'Name', sprintf('%s - duration: %s - averaging: %s', plotSequence.figureName, ...
      get_durationStr(Seq.averageRepetitionTime), get_durationStr(Seq.SequenceTime)));
  else
    if isemptyfield(Seq, 'tOffset')
      tOffset = 0;
    else
      tOffset = Seq.tOffset(1);
    end
    set(hFigure, 'Name', sprintf('%s - duration: %s', plotSequence.figureName, ...
      get_durationStr(tOffset + sum(Seq.tRep(:)))));
  end
  isOwnFigure = true;
elseif ishghandle(hParent, 'uipanel')
  hFigure = ancestor(hParent, 'figure');
  isOwnFigure = false;
else
  error('PD:plotSeq', '"Seq.plotSequence.hParent" must be a valid handle to a figure or uipanel.');
end


%% collect data
% if Seq.plotSeqStart > 1
%   tRepStart = sum(Seq.tRep(1:Seq.plotSeqStart-1));
% else
  tRepStart = 0;
% end
Seq.tRep = Seq.tRep(Seq.plotSeqStart:Seq.plotSeqEnd);
tRepCumsumWraps = reshape(cumsum(Seq.tRep,2), [], plotSequence.wraps);
tStartWraps = [0,tRepCumsumWraps(end,1:end-1)];
if isemptyfield(Seq, 'tOffset')
  warning('PD:plotSeq:noTOffset', ...
    'Seq.tOffset not defined. Assuming 0.');
  Seq.Offset = zeros(1, Seq.plotSeqEnd-Seq.plotSeqStart+1);
else
  Seq.Offset = Seq.tOffset(Seq.plotSeqStart:Seq.plotSeqEnd);
end
OffsetWraps = reshape(Seq.Offset, [], plotSequence.wraps);

if isemptyfield(TX(1), 'Device'), [TX(:).Device] = deal(1); end
hasTX_AmplitudeDC = mod(HW.MMRT(plotSequence.iDevice).FPGA_Firmware, 1e7) >= 0250312;
isPlotTX = false(1, numel(TX)*(hasTX_AmplitudeDC+1));
isPlotTX(plotSequence.plotTX) = true;  % only switch on selected TX structures
isPlotTX([TX(:).Device] ~= plotSequence.iDevice) = false;  % switch off plotting TX for other devices
if hasTX_AmplitudeDC
  % same for DC amplitude in elements > numel(TX)
  isPlotTX(numel(TX)+plotSequence.plotTX) = true;
  isPlotTX([false(1, numel(TX)), [TX(:).Device] ~= plotSequence.iDevice]) = false;
end
TX = TX([TX(:).Device] == plotSequence.iDevice);

talker = PD.Talker.GetInstance();
if ((numel(talker) < plotSequence.iDevice) || talker(plotSequence.iDevice).Dummy...
    || ((HW.MMRT(plotSequence.iDevice).FPGA_PPGeneratorType == 1 ...
         && isprop(talker(plotSequence.iDevice).mySequency.mySequency.myHFPulsesCh0, 'AddChannel'))) ...
        || HW.MMRT(plotSequence.iDevice).FPGA_PPGeneratorType == 2) ...
    && mod(HW.MMRT(plotSequence.iDevice).FPGA_Firmware, 1e7) >= 0210511
  % configuration allows mixing different signals to one channel
  mixedTX = true;
else
  mixedTX = false;
end

for t = 1:length(TX)
  TX(t).Start = TX(t).Start(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  TX(t).Duration = TX(t).Duration(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  TX(t).Amplitude = TX(t).Amplitude(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  TX(t).Frequency = TX(t).Frequency(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  TX(t).Phase = TX(t).Phase(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  if mixedTX && ~isscalar(TX(t).Channel)
    TX(t).Channel = TX(t).Channel(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  end
  if hasTX_AmplitudeDC && ~isemptyfield(TX(t), 'AmplitudeDC')
    TX(t).AmplitudeDC = TX(t).AmplitudeDC(1,Seq.plotSeqStart:Seq.plotSeqEnd);
  end
end

% remove channels with all NaN start from plot
isAllNanTX = cellfun(@(x) all(isnan(x(:))), {TX(:).Start});
TX(isAllNanTX) = [];

if isemptyfield(AQ(1), 'Device'), [AQ(:).Device] = deal(1); end
AQ = AQ([AQ(:).Device] == plotSequence.iDevice);
for t = 1:length(AQ)
  AQ(t).Start = AQ(t).Start(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  AQ(t).nSamples = AQ(t).nSamples(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  if ~isemptyfield(AQ(t), 'SamplingFactor')
    AQ(t).SamplingFactor = AQ(t).SamplingFactor(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  end
  AQ(t).fSample = AQ(t).fSample(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  AQ(t).Frequency = AQ(t).Frequency(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  AQ(t).Phase = AQ(t).Phase(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  if ~isemptyfield(AQ(t), 'FrequencyX')
    AQ(t).FrequencyX = AQ(t).FrequencyX(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  end
  if ~isemptyfield(AQ(t), 'PhaseX')
    AQ(t).PhaseX = AQ(t).PhaseX(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  end
end

if numel(Seq.DigitalIO) >= plotSequence.iDevice
  DigitalIO.SetTime = Seq.DigitalIO(plotSequence.iDevice).SetTime(:,Seq.plotSeqStart:Seq.plotSeqEnd);
  DigitalIO.SetValue = Seq.DigitalIO(plotSequence.iDevice).SetValue(:,Seq.plotSeqStart:Seq.plotSeqEnd);
else
  DigitalIO.SetTime = NaN;
end

PlotTX = double(~isempty(TX) && ~all(isnan(TX(1).Start(:)))) * find(isPlotTX);
PlotDigiIO = plotSequence.plotDigitalIO && any(~isnan(DigitalIO.SetTime(:)));
PlotAQ = ~isempty(AQ) && ~all(isnan(AQ.Start(:)));


%% update figure title
wrapDuration = max(tRepCumsumWraps(end,:)-tStartWraps);
if isOwnFigure && plotSequence.wraps > 1
  set(hFigure, 'Name', sprintf('%s - wrap: %s', get(hFigure, 'Name'), ...
    get_durationStr(wrapDuration)));
end


%% find axes
hAxes = getappdata(hParent, 'plotSeqAxes');
nGrads = HW.Grad(plotSequence.iDevice).n;
if isempty(hAxes), hAxes = zeros(1, nGrads+2); end
hAxesDigiOut = getappdata(hParent, 'plotSeqAxesDigiOut');
hYLabel = getappdata(hParent, 'plotSeqYLabel');
oldPlotSequence = getappdata(hParent, 'plotSequence');
newStart = false;
if isempty(oldPlotSequence)
  oldPlotSequence = plotSequence;
  newStart = true;
end
setappdata(hParent, 'plotSequence', plotSequence);

PlotGrad = false(1, nGrads);
plotSequence.Gradients(plotSequence.Gradients>numel(Grad)) = [];
if plotSequence.Gradients ~= 0
  PlotGrad(plotSequence.Gradients) = true;
end
requiredAxes = [PlotGrad, (any(PlotTX) || PlotDigiIO), PlotAQ];


%% clear foreign elements from parent
hKids = get(hParent, 'Children');
hForeign = findobj(hKids, 'flat', '-not', {'Tag', 'plotSeq_Axes', '-or', 'Tag', 'plotSeqLegend'});
delete(hForeign);


%% prepare axes
if plotSequence.stackGrads && any(PlotGrad)
  if plotSequence.stackTXRX
    numAxes = 1 + double(any(PlotTX) || PlotDigiIO || PlotAQ);
    axesPosition = [PlotGrad, 2, 2];
  else
    numAxes = 1 + sum([(any(PlotTX) || PlotDigiIO), PlotAQ]);
    axesPosition = [PlotGrad, cumsum([(any(PlotTX) || PlotDigiIO), PlotAQ])+1];
  end
elseif plotSequence.stackTXRX
  numAxes = sum(requiredAxes(1:nGrads)) + double(any(PlotTX) || PlotDigiIO || PlotAQ);
  axesPosition = cumsum([requiredAxes(1:nGrads), (any(PlotTX) || PlotDigiIO || PlotAQ)]);
  axesPosition(end+1) = axesPosition(end);
else
  numAxes = sum(requiredAxes);
  axesPosition = cumsum(requiredAxes);
end
doRearrange = false;
if oldPlotSequence.stackGrads ~= plotSequence.stackGrads
  hAxesGrad = hAxes(ishghandle(hAxes(1:nGrads), 'axes'));
  if ~isempty(hAxesGrad)
    if plotSequence.stackGrads
      % find first existing gradient axes and delete all others
      notStackedAxes = hAxes ~= hAxesGrad(1) & [ones(1, nGrads), 0, 0] & ishghandle(hAxes, 'axes');
      delete(hAxes(notStackedAxes));
      hAxes(notStackedAxes) = 0;
      hAxes(PlotGrad) = hAxesGrad(1);
    else
      % re-use "stacked" axes for first requested gradient
      hAxes(1:nGrads) = 0;
      % Always keep one of the references to the (stacked) axes to be able to
      % delete them further down.
      hAxes(max([1, find(PlotGrad, 1, 'first')])) = hAxesGrad(1);
      legend(hAxesGrad(1), 'off');
    end
  end
  doRearrange = true;
end
if oldPlotSequence.stackTXRX ~= plotSequence.stackTXRX
  hAxesTXRX = hAxes([false(1, nGrads), ishghandle(hAxes(5:6), 'axes')]);
  if ~isempty(hAxesTXRX)
    if plotSequence.stackTXRX
      % find first existing TX or RX axes and delete the other
      notStackedAxes = hAxes ~= hAxesTXRX(1) & [zeros(1, nGrads) 1 1] & ishghandle(hAxes, 'axes');
      delete(hAxes(notStackedAxes));
      hAxes(notStackedAxes) = 0;
      hAxes(requiredAxes & [zeros(1, nGrads) 1 1]) = hAxesTXRX(1);
    else
      % re-use "stacked" axes for first requested TX or RX
      hAxes(nGrads+(1:2)) = 0;
      % Always keep one of the references to the (stacked) axes to be able to
      % delete them further down.
      hAxes(max(5, find(requiredAxes & [zeros(1, nGrads) 1 1], 1, 'first'))) = hAxesTXRX(1);
      legend(hAxesTXRX(1), 'off');
    end
  end
  doRearrange = true;
end

if any(ishghandle(hAxes(~requiredAxes), 'axes')) && ...
    (~plotSequence.stackGrads || ~any(PlotGrad))
  % delete axes that are not needed
  delete(hAxes(~requiredAxes & ishghandle(hAxes, 'axes')));
  doRearrange = true;
end
hAxes(~requiredAxes) = 0;

% try to maintain current xlims
prevXLim = [];
hAxesRequired = hAxes(requiredAxes);
presentAxes = ishghandle(hAxesRequired, 'axes');
if any(presentAxes)
  prevXLim = get(hAxesRequired(find(presentAxes, 1, 'first')), 'XLim');
end

% create missing axes
existingAxes = ishghandle(hAxes, 'axes');
reuseAxes = false;
newAxes = false(size(hAxes));
if isempty(hAxes) || any(~ishghandle(hAxes(requiredAxes), 'axes')) || ...
    oldPlotSequence.stackGrads ~= plotSequence.stackGrads || ...
    oldPlotSequence.stackTXRX ~= plotSequence.stackTXRX
  usedStackedGradAxes = false;
  usedStackedTXRXAxes = false;
  for iAxes = 1:numel(requiredAxes)
    decorateAxes = true;
    foundAxes = false;
    if ~requiredAxes(iAxes)
      continue;
    end
    if (((oldPlotSequence.stackGrads == plotSequence.stackGrads || ~usedStackedGradAxes || iAxes > nGrads) || ...
         (oldPlotSequence.stackTXRX  == plotSequence.stackTXRX  || ~usedStackedTXRXAxes || iAxes <= nGrads)) && ...
        existingAxes(iAxes))
      if iAxes <= nGrads && ~usedStackedGradAxes && requiredAxes(iAxes)
        % Re-use the formerly stacked axes for this gradient.
        set(hAxes(iAxes), 'XTickLabel', '');
        xlabel(hAxes(iAxes), '');
        usedStackedGradAxes = true;
      end
      if iAxes > nGrads && ~usedStackedTXRXAxes && requiredAxes(iAxes)
        % Re-use the formerly stacked TX/RX axes.
        set(hAxes(iAxes), 'XTickLabel', '');
        xlabel(hAxes(iAxes), '');
        usedStackedTXRXAxes = true;
      end
      foundAxes = true;
    end
    if ~foundAxes
      decorateAxes = true;
      hGradAxes = hAxes(axesPosition == axesPosition(iAxes));
      % find other axes at this position
      if plotSequence.stackGrads && ...
          any(ishghandle(hGradAxes, 'axes'))
        hAxes(iAxes) = hGradAxes(find(ishghandle(hGradAxes, 'axes'), 1, 'first'));
        foundAxes = true;
        decorateAxes = false;
      end
      if ~foundAxes
        hAxes(iAxes) = axes('Parent', hParent, ...
          'Tag', 'plotSeq_Axes', 'Box', 'on', 'XTickLabel', '');
        newAxes(iAxes) = true;
        doRearrange = true;
      end
    end
    if iAxes <= nGrads
      % decoration for the gradient axes
      t = iAxes+0;
      if decorateAxes
        if plotSequence.stackGrads
          gradLabel = {'Amplitude', ['in ', HW.Grad(plotSequence.iDevice).AmpUnit{t}]};
        else
          gradLabel = {HW.Grad(plotSequence.iDevice).Name{t}, ...
            ['max ', num2str(HW.Grad(plotSequence.iDevice).MaxAmp(t)/HW.Grad(plotSequence.iDevice).AmpUnitScale(t), '%5.3f')], ...
            ['in ', HW.Grad(plotSequence.iDevice).AmpUnit{t}]};
        end
        hYLabel(iAxes) = ylabel(hAxes(iAxes), gradLabel);
      end
    else
      % collect TX output channels used in pulse program
      useTX = true(1, numel(TX));
      for iTX = 1:numel(TX)
        if ~any(TX(iTX).Amplitude(~isnan(TX(iTX).Start))) && ...
            TX(iTX).ChannelOutput ~= HW.TX(plotSequence.iDevice).ChannelDef
          % channel is probably used for (un-)blanking
          useTX(iTX) = false;
        end
      end
      allTXChannels = unique(reshape(cat(1, TX(useTX).Channel), [], 1), 'sorted');
      % FIXME: This does not work for pulses at fLarmorX.
      AmpMax = HW.TX(plotSequence.iDevice).Max.AmplitudeCalibrated(allTXChannels);
      if plotSequence.stackTXRX
        % decoration for the stacked TX/RX axes
        % This label will be overwritten. But it is necessary to set it here to
        % get uniform subplots.
        hYLabel(iAxes) = ylabel(hAxes(iAxes), ...
          {[HW.RX(plotSequence.iDevice).fSampleName, ' in ', HW.RX(plotSequence.iDevice).fSampleUnit], ...
          [HW.TX(plotSequence.iDevice).AmplitudeName ' in ' HW.TX(plotSequence.iDevice).AmplitudeUnit], ...
          sprintf('max %3.1f %s', AmpMax(1)/HW.TX(plotSequence.iDevice).AmplitudeUnitScale, HW.TX(plotSequence.iDevice).AmplitudeUnit)});
      else
        if iAxes == nGrads+1
          % decoration for the TX axes
          % This label will be overwritten. But it is necessary to set it here to
          % get uniform subplots.
          % FIXME: It should be possible to use different amplitude names, units
          % and unit scales for signals on Tx1 and Tx2.
          TXYLabelString = cell(numel(AmpMax)+1, 1);
          TXYLabelString{1} = [HW.TX(plotSequence.iDevice).AmplitudeName, ' in ', HW.TX(plotSequence.iDevice).AmplitudeUnit];
          for iTX = 1:numel(allTXChannels)
            TXYLabelString{iTX+1} = sprintf('max %3.1f %s (ch %d)', ...
              AmpMax(iTX)/HW.TX(plotSequence.iDevice).AmplitudeUnitScale, ...
              HW.TX(plotSequence.iDevice).AmplitudeUnit, ...
              allTXChannels(iTX));
          end
          hYLabel(iAxes) = ylabel(hAxes(nGrads+1), TXYLabelString);
        else % if iAxes == 6
          % decoration for RX axes
          hYLabel(iAxes) = ylabel(hAxes(nGrads+2), {HW.RX(plotSequence.iDevice).fSampleName, ['in ', HW.RX(plotSequence.iDevice).fSampleUnit]});
        end
      end
    end
    % decoration for all axes
    grid(hAxes(iAxes), 'on');
    set(hAxes(iAxes), 'XMinorGrid', 'on', 'YMinorGrid', 'on');
  end

  setappdata(hParent, 'plotSeqAxes', hAxes);
  setappdata(hParent, 'plotSeqYLabel', hYLabel);
else
  % found exact number of axes in parent
  reuseAxes = true;
end

% re-arrange all axes
if doRearrange
  % To get uniform behavior, do this after all axes have "ylabel" set
  % FIXME: This still doesn't work for the lowest axes. But the non-uniformity
  % is only visible when scaling the figure very narrow.
  % linkaxes(hAxes(requiredAxes), 'off');
  for iAxes = 1:numel(requiredAxes)
    if requiredAxes(iAxes)
      % xlabel(hAxes(iAxes), '')
      % set(hAxes(iAxes), 'XTickLabel', '', 'ActivePositionProperty', 'position', ...
      %  'Position', get(hAxes(iAxes), 'Position'), 'XLimMode', 'auto');
      hAxes(iAxes) = subplot(numAxes,1,axesPosition(iAxes), hAxes(iAxes));
    end
  end
  if any(ishghandle(hAxes(requiredAxes)))
    linkaxes(hAxes(requiredAxes), 'x');
    % disable rotate3d for these axes
    hrot3d = rotate3d(hFigure);
    hrot3d.setAllowAxesRotate(hAxes(requiredAxes), false);
  end
end

if plotSequence.stackGrads && any(PlotGrad)
  % add legend to plot
  gradAxes = find(PlotGrad);
  hLineLegend = findobj(get(hAxes(gradAxes(1)), 'Children'), 'flat', 'Tag', 'plotSeqLineLegend');
  for iAxes = gradAxes
    t = iAxes+0;
    % insert "invisible" line with center color for legend
    if numel(hLineLegend) < sum(requiredAxes(1:iAxes))
      hLineLegend = ...
        [line('Parent', hAxes(iAxes), 'XData', NaN, 'YData', NaN, 'Color', mean(GradColors{t}, 1), 'Tag', 'plotSeqLineLegend');
        hLineLegend];
      % sort legend line to front
      hOtherKids = findobj(get(hAxes(iAxes), 'Children'), 'flat', '-not', 'Tag', 'plotSeqLineLegend');
      set(hAxes(iAxes), 'Children', [hOtherKids; hLineLegend]);
    else
      set(hLineLegend(end-sum(requiredAxes(1:iAxes))+1), 'Color', mean(GradColors{t}, 1));
    end
  end
  gradLegend = HW.Grad(plotSequence.iDevice).Name(gradAxes);
  if newStart || (oldPlotSequence.stackGrads ~= plotSequence.stackGrads) || ...
      ~isequal(plotSequence.Gradients, oldPlotSequence.Gradients)
    % decoration for common gradient axes
    hCommonAxes = hAxes(ishghandle(hAxes(1:nGrads), 'axes'));
    lh = legend(hCommonAxes(1), gradLegend, 'Tag', 'plotSeqLegend', ...
      'Orientation', 'horizontal', 'Location', 'north');
    legend(hCommonAxes(1), 'boxoff');
    if ~verLessThan('Matlab', '9.2')
      set(lh, 'AutoUpdate', 'off'); % Do not automatically add new entries to legend
    end
    % move legend outside manually to keep the axes height
    legendPos = get(lh, 'Position');
    axPos = get(hCommonAxes(1), 'Position');
    legendPos(2) = axPos(2) + axPos(4);
    set(lh, 'Position', legendPos);
  end
end

% axes for Digital IO
if PlotDigiIO  % any(~isnan(DigitalIO.SetTime(:)))
  if isempty(hAxesDigiOut) || ~ishghandle(hAxesDigiOut, 'axes')
    hAxesDigiOut = axes('Parent', hParent, 'Position', get(hAxes(nGrads+1), 'Position'), ...
      'YAxisLocation', 'right', 'Color', 'none', 'Box', 'off', 'HitTest', 'off', ...
      'XTickLabel', '', 'XGrid', 'off', 'YGrid', 'off', ...
      'YLim', [0.5, 6.5], 'YTick', 1:6, 'YGrid', 'off', ...
      'YColor', plotSequence.DigitalIOColor, 'Tag', 'plotSeq_Axes', 'TickLength', [0, 0]);
    ylabel(hAxesDigiOut, 'Digital Output');

    hLink = linkprop([hAxes(nGrads+1), hAxesDigiOut], 'Position');
    setappdata(hAxesDigiOut, 'plotSeqAxesDigiOutPositionLinker', hLink);

    linkaxes([hAxes(requiredAxes), hAxesDigiOut], 'x');
    doRearrange = true;
  end
else
  if ishghandle(hAxesDigiOut, 'axes')
    delete(hAxesDigiOut);
  end
  hAxesDigiOut = [];
end
setappdata(hParent, 'plotSeqAxesDigiOut', hAxesDigiOut);

% set font sizes
if ~isempty(plotSequence.FontSize)
  set(hAxes(requiredAxes), 'FontSize', plotSequence.FontSize, ...
    'LabelFontSizeMultiplier', 1);
  if ishghandle(hAxesDigiOut, 'axes')
    set(hAxesDigiOut, 'FontSize', plotSequence.FontSize, ...
    'LabelFontSizeMultiplier', 1);
  end
end

% Update xlims
lastValidAxes = hAxes(find(ishghandle(hAxes, 'axes'), 1, 'last'));
if ~isemptyfield(plotSequence, 'xLim'), xLim = plotSequence.xLim; else xLim = [-Inf, Inf]; end
if isinf(xLim(1)), xLim(1) = -max(OffsetWraps(1,:))+tRepStart; end
if isinf(xLim(2)), xLim(2) = max(tRepCumsumWraps(end,:)-tStartWraps)+tRepStart; end

oldXLim = getappdata(hParent, 'plotSeqXLim');
setappdata(hParent, 'plotSeqXLim', xLim);


%% plot sequence
if plotSequence.Gradients ~= 0
  numLinesStacked = 0;
  maxAmp = 0;
  for t = plotSequence.Gradients
    % plot tRep connections und nan areas
    shimOffset = HW.Grad(plotSequence.iDevice).AmpOffsetExtra(t) + ...
      (Grad(t).Shim + HW.Grad(plotSequence.iDevice).AmpOffset(t)) * plotSequence.includeShim;
    Grad(t).Amp = Grad(t).Amp(:,Seq.plotSeqStart:Seq.plotSeqEnd) + shimOffset;
    Grad(t).Time = Grad(t).Time(:,Seq.plotSeqStart:Seq.plotSeqEnd);

    firstAmp = Grad(t).Amp(1,:);
    startTime = zeros(size(Seq.tRep)) - Seq.Offset ...
      + HW.Grad(plotSequence.iDevice).TimeDelay(HW.Grad(plotSequence.iDevice).xyzB(t));
    lastAmpindex = sum(~isnan(Grad(t).Time), 1);
    lastAmpindex0 = (lastAmpindex==0);
    lastAmpindex(lastAmpindex==0) = 1;
    lastAmp = Grad(t).Amp(lastAmpindex + size(Grad(t).Amp,1)*(0:size(Grad(t).Amp,2)-1));
    lastTime = Grad(t).Time(lastAmpindex + size(Grad(t).Amp,1)*(0:size(Grad(t).Amp,2)-1));
    if plotSequence.wraps==1
      endTime = Seq.tRep ...
        - [Seq.Offset(2:end), (double((Seq.averageBreak==0)&&(Seq.average>1)) * Seq.Offset(1))] ...
        + HW.Grad(plotSequence.iDevice).TimeDelay(HW.Grad(plotSequence.iDevice).xyzB(t));
    else
      endTime = Seq.tRep ...
        - [Seq.Offset(2:end), NaN] ...
        + HW.Grad(plotSequence.iDevice).TimeDelay(HW.Grad(plotSequence.iDevice).xyzB(t));
    end
    lastAmp(lastAmpindex0) = shimOffset;
    lastTime(lastAmpindex0) = startTime(lastAmpindex0);
    firstAmp(lastAmpindex0) = shimOffset;

    tt = ones(size(Grad(t).Time,1)+5,1)*cumsum([0,Seq.tRep(1:end-1)]);
    % plot([[-Seq.Offset(1),startTime(1:end-1)]; startTime;            startTime; Grad(t).Time;lastTime;endTime]+tt+tRepStart,...
    %      [[0,lastAmp(1:end-1)];                [0,lastAmp(1:end-1)]; firstAmp;  Grad(t).Amp; lastAmp; lastAmp]./HW.Grad(plotSequence.iDevice).AmpUnitScale(t))
    % plot(hax(at), ...
    %   [startTime;            startTime;            startTime; Grad(t).Time; lastTime; endTime]+tt+tRepStart, ...
    %   [[0,lastAmp(1:end-1)]; [0,lastAmp(1:end-1)]; firstAmp;  Grad(t).Amp;  lastAmp;  lastAmp]./HW.Grad(plotSequence.iDevice).AmpUnitScale(t));

    wrappedTime = round(bsxfun(@minus, reshape([startTime;                      startTime;                    startTime; Grad(t).Time; lastTime; endTime]+tt+tRepStart,             [], plotSequence.wraps), ...
                               tStartWraps)*HW.MMRT(plotSequence.iDevice).fSystem)/HW.MMRT(plotSequence.iDevice).fSystem;
    wrappedAmp  =                      reshape([[shimOffset,lastAmp(1:end-1)]; [shimOffset,lastAmp(1:end-1)]; firstAmp;  Grad(t).Amp;  lastAmp;  lastAmp]./HW.Grad(plotSequence.iDevice).AmpUnitScale(t), [], plotSequence.wraps);

    uniqueTime = wrappedTime; uniqueTime(isnan(uniqueTime)) = 1/eps;
    uniqueAmp = wrappedAmp;   uniqueAmp(isnan(uniqueAmp)) = 1/eps;

    if all(all(abs(bsxfun(@minus, uniqueTime(1:end,1), uniqueTime(1:end,:))) < 3*eps)) && ...
        all(all(abs(bsxfun(@minus, uniqueAmp(1:end,1), uniqueAmp(1:end,:))) < 3*eps))
      colorOrder = {mean(GradColors{t}, 1)};
      uniqueTime = wrappedTime(:,1);
      uniqueAmp = wrappedAmp(:,1);
    else
      roundTime = wrappedTime - min(wrappedTime(:));
      roundTime = int32(roundTime/wrapDuration*(2^31-1));
      roundAmp = int32(wrappedAmp/max(HW.Grad(plotSequence.iDevice).MaxAmp(t), max(abs(wrappedAmp(:))))*(2^31-1));
      [~, iA] = unique([roundTime; roundAmp].', 'rows', 'stable');

      uniqueTime = wrappedTime(:,iA);
      uniqueAmp = wrappedAmp(:,iA);
      numColors = size(GradColors{t}, 1);
      colorOrder = mat2cell(interp1(linspace(0, 1, numColors), GradColors{t}, ...
        linspace(0, 1, size(uniqueTime, 2))), ones(1, size(uniqueTime, 2)), 3);
    end

    uniqueTime = uniqueTime - plotSequence.tOffset;

    if ~plotSequence.stackGrads || t == find(PlotGrad, 1, 'first')
      % find lines in axes
      hKids = get(hAxes(t), 'Children');
      hLines = flipud(findobj(hKids, 'flat', 'Type', 'line', '-not', 'Tag', 'plotSeqLineLegend'));
    end
    hLinesCurr = hLines;
    % compare number of found lines with required lines
    numLinesReq = size(uniqueAmp, 2);
    numLinesPre = numel(hLines) - numLinesStacked;
    if numLinesPre > numLinesReq
      if ~plotSequence.stackGrads || t == find(PlotGrad, 1, 'last')
        % delete surplus lines
        delete(hLines(numLinesStacked+numLinesReq+1:end));
      end
      hLinesCurr(numLinesStacked+numLinesReq+1:end) = [];
      numLinesPre = numLinesReq;
    end
    uniqueTimePre = mat2cell(uniqueTime(:,1:numLinesPre).', ones(1, numLinesPre), size(uniqueTime, 1));
    uniqueAmpPre = mat2cell(uniqueAmp(:,1:numLinesPre).', ones(1, numLinesPre), size(uniqueAmp, 1));
    set(hLinesCurr(numLinesStacked+1:end), {'XData'}, uniqueTimePre, {'YData'}, uniqueAmpPre, 'Tag', sprintf('plotSeqLineGrad%d', t));
    % add missing lines
    hLinesNew = line(uniqueTime(:,numLinesPre+1:end), uniqueAmp(:,numLinesPre+1:end), 'Tag', sprintf('plotSeqLineGrad%d', t), 'Parent', hAxes(t));
    hLinesCurr = [hLinesCurr(numLinesStacked+1:end); hLinesNew];
    hLines = [hLines; hLinesNew];
    set(hLinesCurr, {'Color'}, colorOrder);
    if plotSequence.stackGrads
      numLinesStacked = numLinesStacked + numel(hLinesCurr);
      maxAmp = max(maxAmp, max(abs(uniqueAmp(~isnan(uniqueAmp)))));
    else
      maxAmp = max(abs(uniqueAmp(~isnan(uniqueAmp))));
      if maxAmp > 0
        set(hAxes(t), 'YLim', [-maxAmp, maxAmp]*1.1);
      end
    end
  end
  if plotSequence.stackGrads && maxAmp > 0
    yLim = [-maxAmp, maxAmp]*1.1;
    oldYLim = getappdata(hParent, 'plotSeqGradYLim');
    setappdata(hParent, 'plotSeqGradYLim', yLim);
    if ~isequal(oldYLim, yLim)
      set(hAxes(PlotGrad), 'YLim', yLim);
    end
  else
    setappdata(hParent, 'plotSeqGradYLim', [0 0]);
  end
end

numLinesTX = 0;
if any(PlotTX)
  tempmax = 0;
  tempmin = 0;
  oldTXLength = getappdata(hParent, 'plotSeqTXLength');
  if isempty(oldTXLength), oldTXLength = 0; end
  TXLength = numel(PlotTX);
  setappdata(hParent, 'plotSeqTXLength', TXLength);
  legendCell = cell(1, TXLength);

  % find lines in axes
  hKids = get(hAxes(nGrads+1), 'Children');
  hLines = flipud(findobj(hKids, 'flat', 'Type', 'line', '-not', 'Tag', 'plotSeqLineLegend'));
  legendHandles = gobjects(1, TXLength);

  for iTX = reshape(PlotTX, 1, [])
    if iTX > numel(TX)
      tt = iTX - numel(TX);
      useAmplitudeDC = true;
    else
      tt = iTX;
      useAmplitudeDC = false;
    end
    if ~any(TX(tt).Amplitude(~isnan(TX(tt).Start))) && ...
        TX(tt).ChannelOutput ~= HW.TX(plotSequence.iDevice).ChannelDef
      % channel is probably used for (un-)blanking
      TXLength = TXLength - 1;
      continue;
    end

    % FIXME: This doesn't take the actual output channel into account

    %  arrays for plotting
    t = ones(size(TX(tt).Start,1),1) * cumsum([0,Seq.tRep(1:end-1)]);
    TX(tt).tTxPlot = reshape([t(:).'+TX(tt).Start(:).'; ...
                              t(:).'+TX(tt).Start(:).'; ...
                              t(:).'+TX(tt).Start(:).'+TX(tt).Duration(:).'; ...
                              t(:).'+TX(tt).Start(:).'+TX(tt).Duration(:).'; ...
                              ], [size(TX(tt).Start,1)*4, size(TX(tt).Start,2), size(TX(tt).Start,3)])+tRepStart;
    if useAmplitudeDC
      isTXidx = isfinite(TX(tt).Amplitude);
      TX(tt).Amplitude = repmat(TX(tt).AmplitudeDC, size(TX(tt).Amplitude, 1), 1);
      TX(tt).Phase(:) = 0;
      TX(tt).Frequency(:) = 0;
      TX(tt).Amplitude(~isTXidx) = NaN;
      if all(TX(tt).Amplitude(isTXidx) == 0)
        % all DC amplitudes are zero, so don't plot a line for it
        TXLength = TXLength - 1;
        continue;
      end
    end
    TX(tt).AmplitudeTxPlot = (1/HW.TX(plotSequence.iDevice).AmplitudeUnitScale)*reshape([ ...
                                  inf(size(TX(tt).Amplitude(:),1),1).'; ...
                                  TX(tt).Amplitude(:).'; ...
                                  TX(tt).Amplitude(:).'; ...
                                  inf(size(TX(tt).Amplitude(:),1),1).' ...
                                  ], [size(TX(tt).Amplitude,1)*4, size(TX(tt).Amplitude,2), size(TX(tt).Amplitude,3)]);
    TX(tt).PhaseTxPlot = reshape([zeros(size(TX(tt).Phase(:),1),1).'; ...
                                  TX(tt).Phase(:).'; ...
                                  TX(tt).Phase(:).'; ...
                                  zeros(size(TX(tt).Phase(:),1),1).' ...
                                  ], [size(TX(tt).Phase,1)*4, size(TX(tt).Phase,2), size(TX(tt).Phase,3)]);

    TX(tt).AmpComplexTxPlot = TX(tt).AmplitudeTxPlot.*(cosd(TX(tt).PhaseTxPlot) + 1i*sind(TX(tt).PhaseTxPlot));

    TX(tt).FrequencyTxPlot = reshape([zeros(size(TX(tt).Frequency(:),1),1).'; ...
                                  TX(tt).Frequency(:).'; ...
                                  TX(tt).Frequency(:).'; ...
                                  zeros(size(TX(tt).Frequency(:),1),1).' ...
                                  ], [size(TX(tt).Frequency,1)*4, size(TX(tt).Frequency,2), size(TX(tt).Frequency,3)]);

    if mixedTX
      TX(tt).ChannelTxPlot = reshape(repmat(TX(tt).Channel(:).', 4, 1), ...
        [size(TX(tt).Channel,1)*4, size(TX(tt).Channel,2), size(TX(tt).Channel,3)]);
    end

    % wrap
    ampTX = reshape(TX(tt).AmpComplexTxPlot(:), [], plotSequence.wraps);
    fSystem = HW.MMRT(plotSequence.iDevice).fSystem;
    tTX = round(bsxfun(@minus, reshape(round(TX(tt).tTxPlot(:)*fSystem)/fSystem, [], plotSequence.wraps), tStartWraps)*fSystem)/fSystem;
    freqTX = reshape(TX(tt).FrequencyTxPlot(:), [], plotSequence.wraps);
    if mixedTX
      channelTX = reshape(TX(tt).ChannelTxPlot(:), [], plotSequence.wraps);
    end

    % remove zeros at the same time
    sameTime = abs(tTX(1:end-1,:)-tTX(2:end,:)) < 1.5/fSystem;
    ampInf = isinf(ampTX);
    toRemove = ampInf & ... % has zero amplitude
               [false(1,size(tTX,2)); sameTime] & ... % and same time as before
               [sameTime; false(1,size(tTX,2))]; % and same time as after
    tempIndex = bsxfun(@plus, cumsum(~toRemove, 1), (0:(size(toRemove,2)-1))*size(toRemove,1));
    tempIndex(toRemove) = NaN;
    ampTX1 = NaN(size(tTX));
    tTX1 = NaN(size(tTX));
    freqTX1 = NaN(size(tTX));
    ampTX1(tempIndex(~isnan(tempIndex))) = ampTX(~toRemove);
    tTX1(tempIndex(~isnan(tempIndex))) = tTX(~toRemove);
    freqTX1(tempIndex(~isnan(tempIndex))) = freqTX(~toRemove);
    nansRemove = all(isnan(tTX1), 2);
    ampTX1(nansRemove,:) = [];
    tTX1(nansRemove,:) = [];
    freqTX1(nansRemove,:) = [];
    if mixedTX
      channelTX1 = NaN(size(tTX));
      channelTX1(tempIndex(~isnan(tempIndex))) = channelTX(~toRemove);
      channelTX1(nansRemove,:) = [];
    end

    % insert nans between Inf
    ampInf = isinf(ampTX1);
    insertAt = ampInf & [false(1,size(ampTX1,2)); ampInf(1:end-1,:)];
    numRowsAdd = max(sum(insertAt,1));
    tempIndex = bsxfun(@plus, bsxfun(@plus, cumsum(insertAt,1), (1:size(insertAt,1)).'), ...
      (0:(size(insertAt,2)-1))*(size(insertAt,1)+numRowsAdd));
    tTX = nan(size(tTX1,1)+numRowsAdd, size(tTX1,2));
    ampTX = tTX;
    freqTX = tTX;
    if mixedTX
      channelTX = tTX;
    end
    tTX(tempIndex(:)) = tTX1(:);
    ampTX(tempIndex(:)) = ampTX1(:);
    freqTX(tempIndex(:)) = freqTX1(:);
    if mixedTX
      channelTX(tempIndex(:)) = channelTX1(:);
    end
    ampTX(isinf(ampTX))=0; % replace inf with zeros

    % check if all wraps are equal
    uniqueTime = tTX; uniqueTime(isnan(uniqueTime)) = 1/eps;
    uniqueAmp = ampTX; uniqueAmp(isnan(uniqueAmp)) = 1/eps;
    uniqueFreq = freqTX; uniqueFreq(isnan(uniqueFreq)) = 1/eps;
    if mixedTX
      uniqueChannel = channelTX; uniqueChannel(isnan(uniqueChannel)) = 1/eps;
    end
    if all(all(abs(bsxfun(@minus, uniqueTime(:,1), uniqueTime)) < 3*eps)) && ...
        all(all(abs(bsxfun(@minus, uniqueAmp(:,1), uniqueAmp)) < 3*eps)) && ...
        all(all(abs(bsxfun(@minus, uniqueFreq(:,1), uniqueFreq)) < 3*eps)) && ...
        (mixedTX &&all(all(abs(bsxfun(@minus, uniqueChannel(:,1), uniqueChannel)) < 3*eps)))
      uniqueTime = tTX(:,1);
      uniqueAmp = ampTX(:,1);
      uniqueFreq = freqTX(:,1);
      if mixedTX
        uniqueChannel = channelTX(:,1);
      end
    else
      roundTime = tTX - min(tTX(:));
      roundTime = int32(roundTime/wrapDuration*(2^31-1));
      roundRealAmp = int32(real(ampTX)/HW.TX(plotSequence.iDevice).AmpMax*(2^31-1));
      roundImagAmp = int32(imag(ampTX)/HW.TX(plotSequence.iDevice).AmpMax*(2^31-1));
      [~, iA] = unique([roundTime; roundRealAmp; roundImagAmp].', 'rows', 'stable');

      uniqueTime = tTX(:,iA);
      uniqueAmp = ampTX(:,iA);
      uniqueFreq = freqTX(:,iA);
      if mixedTX
        uniqueChannel = channelTX(:,iA);
      end
    end

    uniqueTime = uniqueTime - plotSequence.tOffset;

    % compare number of found lines with required lines
    numLinesReq = size(uniqueAmp, 2) * double(size(uniqueAmp, 1) > 0);
    numLinesPre = max(0, numel(hLines)-numLinesTX);
    uniqueTimePre = mat2cell(uniqueTime, size(uniqueTime, 1), ones(1, size(uniqueTime, 2))).';
    uniqueAmpPre = mat2cell(uniqueAmp, size(uniqueAmp, 1), ones(1, size(uniqueAmp, 2))).';
    uniqueFreqPre = mat2cell(uniqueFreq, size(uniqueFreq, 1), ones(1, size(uniqueFreq, 2))).';
    if mixedTX
      uniqueChannelPre = mat2cell(uniqueChannel, size(uniqueChannel, 1), ones(1, size(uniqueChannel, 2))).';
      uniqueUserData = cellfun(@(x,y,z) struct('abs', abs(x), 'phase', angle(x), 'freq', y, 'channel', z), ...
        uniqueAmpPre, uniqueFreqPre, uniqueChannelPre, 'UniformOutput', false);
    else
      uniqueUserData = cellfun(@(x,y) struct('abs', abs(x), 'phase', angle(x), 'freq', y), ...
        uniqueAmpPre, uniqueFreqPre, 'UniformOutput', false);
    end

    % replace or add new lines
    if useAmplitudeDC
      tagTX_str = 'DC';
      lineStyle = ':';
      lineWidth = 2.5;
      absFcn = @(x) x;  % no-op, i.e., allow negative values without "phase"
    else
      tagTX_str = '';
      lineStyle = '-';
      lineWidth = 2;
      absFcn = @abs;
    end
    if mixedTX
      tagTX = sprintf('plotSeqLineTX%d%s', tt, tagTX_str);
    else
      tagTX = sprintf('plotSeqLineTX%d%s', TX(tt).Channel, tagTX_str);
    end
    lastLine = min(numLinesReq, numLinesPre);
    set(hLines(numLinesTX+(1:lastLine)), {'XData'}, uniqueTimePre(1:lastLine), ...
      {'YData'}, cellfun(absFcn, uniqueAmpPre(1:lastLine), 'UniformOutput', false), ...
      'ZData', [], {'UserData'}, uniqueUserData(1:lastLine), ...
      'Color', plotSequence.TXColors{tt,1}, 'LineWidth', lineWidth, 'LineStyle', lineStyle, ...
      'Tag', tagTX);
    hl = line(uniqueTime(:,lastLine+1:end), ...
      absFcn(uniqueAmp(:,lastLine+1:end)), ...
      'Color', plotSequence.TXColors{tt,1}, 'LineWidth', lineWidth, 'LineStyle', lineStyle, ...
      'Tag', tagTX, 'Parent', hAxes(nGrads+1));
    set(hl, {'UserData'}, uniqueUserData(lastLine+1:end));
    currentLineHandles = [hLines(numLinesTX+(1:lastLine)); hl];
    legendHandles(iTX) = currentLineHandles(1);
    if plotSequence.showTXPhase && ~useAmplitudeDC
      lastLine2 = min(2*numLinesReq, numLinesPre);
      set(hLines(numLinesTX+(lastLine+1:lastLine2)), {'XData'}, uniqueTimePre(1:(lastLine2-numLinesReq)), ...
        {'YData'}, cellfun(@real, uniqueAmpPre(1:(lastLine2-numLinesReq)), 'UniformOutput', false), ...
        'ZData', [], {'UserData'}, uniqueUserData(1:(lastLine2-numLinesReq)), ...
        'Color', plotSequence.TXColors{tt,2}, 'LineWidth', 0.5, 'LineStyle', '-', ...
        'Tag', tagTX);
      hl = line(uniqueTime(:,max(1,lastLine2-numLinesReq+1):end), ...
        real(uniqueAmp(:,max(1,lastLine2-numLinesReq+1):end)), ...
        'Color', plotSequence.TXColors{tt,2}, 'LineWidth', 0.5, 'LineStyle', '-', ...
        'Tag', tagTX, 'Parent', hAxes(nGrads+1));
      set(hl, {'UserData'}, uniqueUserData(max(1,lastLine2-numLinesReq+1):end));
      lastLine3 = min(3*numLinesReq, numLinesPre);
      set(hLines(numLinesTX+(lastLine2+1:lastLine3)), {'XData'}, uniqueTimePre(1:(lastLine3-2*numLinesReq)), ...
        {'YData'}, cellfun(@imag, uniqueAmpPre(1:(lastLine3-2*numLinesReq)), 'UniformOutput', false), ...
        'ZData', [], {'UserData'}, uniqueUserData(1:(lastLine3-2*numLinesReq)), ...
        'Color', plotSequence.TXColors{tt,3}, 'LineWidth', 0.5, 'LineStyle', '-', ...
        'Tag', tagTX);
      hl = line(uniqueTime(:,max(1,lastLine3-2*numLinesReq+1):end), ...
        imag(uniqueAmp(:,max(1,lastLine3-2*numLinesReq+1):end)), ...
        'Color', plotSequence.TXColors{tt,3}, 'LineWidth', 0.5, 'LineStyle', '-', ...
        'Tag', tagTX, 'Parent', hAxes(nGrads+1));
      set(hl, {'UserData'}, uniqueUserData(max(1,lastLine3-2*numLinesReq+1):end));
    end

    if useAmplitudeDC
      tagTX_str = ' DC';
    else
      tagTX_str = '';
    end
    legendCell{iTX} = ['TX Channel ', num2str(TX(tt).ChannelOutput), tagTX_str];
    tempmax = max(tempmax, max(abs(TX(tt).AmpComplexTxPlot(:))));
    tempmin = min(tempmin, min(  min(imag(TX(tt).AmpComplexTxPlot(:))),    min(real(TX(tt).AmpComplexTxPlot(:)))  ));
    numLinesTX = numLinesTX + (1+2*(plotSequence.showTXPhase && ~useAmplitudeDC))*numLinesReq;
  end

  % Update TXLength in graphics object. Drawing some lines might have been
  % skipped in the above loop.
  setappdata(hParent, 'plotSeqTXLength', TXLength);

  if numel(hLines) > numLinesTX
    if ~plotSequence.stackTXRX
      % delete surplus lines
      delete(hLines(numLinesTX+1:end));
    end
  end

  if TXLength ~= oldTXLength
    if TXLength > 1 || ...
        any(TX(1).Channel(~isnan(TX(1).Channel)) ~= HW.TX(plotSequence.iDevice).ChannelDef)
      % Avoid issue if the data in TX(1) was identified as a blanking signal and
      % no line was plotted for it.
      validHandles = isgraphics(legendHandles);
      legend(legendHandles(validHandles), legendCell(validHandles), 'Tag', 'plotSeqLegend');
      legend(hAxes(nGrads+1), 'boxoff');
    else
      legend(hAxes(nGrads+1), 'off');
    end
  end

  % Re-set the axes decorations. They might have been deleted if a pulse program
  % included digital outputs but no rf pulses.
  % collect TX output channels used in pulse program
  useTX = true(1, numel(TX));
  for iTX = 1:numel(TX)
    if ~any(TX(iTX).Amplitude(~isnan(TX(iTX).Start))) && ...
        TX(iTX).ChannelOutput ~= HW.TX(plotSequence.iDevice).ChannelDef
      % channel is probably used for (un-)blanking
      useTX(iTX) = false;
    end
  end
  allTXChannels = unique(reshape(cat(1, TX(useTX).Channel), [], 1), 'sorted');
  % FIXME: This does not work for pulses at fLarmorX.
  AmpMax = HW.TX(plotSequence.iDevice).Max.AmplitudeCalibrated(allTXChannels);
  if plotSequence.stackTXRX
    TXYLabelString = {[HW.RX(plotSequence.iDevice).fSampleName, ' in ', HW.RX(plotSequence.iDevice).fSampleUnit], ...
      [HW.TX(plotSequence.iDevice).AmplitudeName ' in ' HW.TX(plotSequence.iDevice).AmplitudeUnit], ...
      sprintf('max(%d) %3.1f %s', allTXChannels(1), AmpMax(1)/HW.TX(plotSequence.iDevice).AmplitudeUnitScale, HW.TX(plotSequence.iDevice).AmplitudeUnit)};
  else
    TXYLabelString = cell(numel(AmpMax)+1, 1);
    TXYLabelString{1} = [HW.TX(plotSequence.iDevice).AmplitudeName, ' in ', HW.TX(plotSequence.iDevice).AmplitudeUnit];
    for iTX = 1:numel(allTXChannels)
      TXYLabelString{iTX+1} = sprintf('max %3.1f %s (ch %d)', ...
        AmpMax(iTX)/HW.TX(plotSequence.iDevice).AmplitudeUnitScale, ...
        HW.TX(plotSequence.iDevice).AmplitudeUnit, ...
        allTXChannels(iTX));
    end
  end
  oldTXYLabelString = getappdata(hAxes(nGrads+1), 'plotSeqYLabel');
  setappdata(hAxes(nGrads+1), 'plotSeqYLabel', TXYLabelString);
  if ~isequal(TXYLabelString, oldTXYLabelString)
    set(get(hAxes(nGrads+1), 'YLabel'), 'String', TXYLabelString);
  end

  if plotSequence.showTXPhase
    if (tempmin ~= tempmax)
      set(hAxes(nGrads+1), 'YLim', [tempmin*1.1, tempmax*1.1]);
    end
  elseif tempmax ~= 0
    set(hAxes(nGrads+1), 'YLim', [0, tempmax*1.1]);
  end
  set(hAxes(nGrads+1), 'YTickMode', 'auto');

elseif PlotDigiIO && ishghandle(hAxes(nGrads+1))
  delete(get(hAxes(nGrads+1), 'Children'));
end

if PlotDigiIO
  if ~any(PlotTX)
    set(hYLabel(nGrads+1), 'String', '');
    set(hAxes(nGrads+1), 'YTick', []);
  end
  % Plot Digital IO
  tIO = DigitalIO.SetTime + repmat(cumsum([0,Seq.tRep(1:end-1)]),size(DigitalIO.SetTime,1),1) + tRepStart;
  tIO = round(bsxfun(@minus, reshape(tIO, [], plotSequence.wraps), tStartWraps)*HW.MMRT(plotSequence.iDevice).fSystem)/HW.MMRT(plotSequence.iDevice).fSystem;
  digiOut = reshape(DigitalIO.SetValue, [], plotSequence.wraps);

  % insert point at start with state from previous wrap and at stop with last
  tIO     = [NaN(1,size(tIO,2));     tIO;     NaN(1,size(tIO,2))];
  digiOut = [NaN(1,size(digiOut,2)); digiOut; NaN(1,size(digiOut,2))];
  tIO(1,:) = -OffsetWraps(1,:);
  if plotSequence.wraps == 1
    tIO(end) = tRepCumsumWraps(end) - double((Seq.averageBreak==0) && (Seq.average>1)) * OffsetWraps(1);
  else
    tIO(end,1:end-1) = diff(tStartWraps) - OffsetWraps(1,2:end);
  end
  digiOut(end,:) = digiOut(end-1,:); digiOut(1,2:end) = digiOut(end,1:end-1);

  % check if all wraps are equal
  uniqueTime = tIO; uniqueTime(isnan(uniqueTime)) = 1/eps;
  uniqueDigiOut = digiOut; uniqueDigiOut(isnan(uniqueDigiOut)) = double(intmax('int32'));
  if all(all(abs(bsxfun(@minus, uniqueTime(:,1), uniqueTime)) < 3*eps)) && ...
      all(all(abs(bsxfun(@minus, uniqueDigiOut(:,1), uniqueDigiOut)) < 3*eps))
    uniqueTime = tIO(:,1);
    uniqueDigiOut = digiOut(:,1);
  else
    roundTime = tIO - min(tIO(:));
    roundTime = int32(roundTime/wrapDuration*(2^31-1));
    [~, iA] = unique([roundTime; int32(uniqueDigiOut)].', 'rows', 'stable');

    uniqueTime = tIO(:,iA);
    uniqueDigiOut = digiOut(:,iA);
  end

  % remove or fill in NaNs
  uniqueDigiOut(all(isnan(uniqueTime), 2),:) = [];
  uniqueTime(all(isnan(uniqueTime), 2),:) = [];
  uniqueDigiOut = fillNaNs(uniqueDigiOut, 1);
  uniqueTime = fillNaNs(uniqueTime, 1);

  uniqueTime = uniqueTime - plotSequence.tOffset;

  % separate channels
  uniqueDigiOut(isnan(uniqueDigiOut)) = 0; % FIXME: Is there a better way to handle NaNs?
  channels = log2(bsxfun(@bitand, uniqueDigiOut, reshape(2.^(0:5), 1, 1, [])))+1;
  channels(isinf(channels)) = 0;
  % insert point after end of signal high (per channel)
  diffChannels = cat(1, diff(channels), zeros(1, size(channels, 2), size(channels, 3)));
  insertRow = any(any(diffChannels<0, 2), 3);
  channels2 = NaN(size(channels, 1) + sum(insertRow), size(channels, 2), size(channels, 3));
  tmp = cumsum(insertRow);
  insertRowAt = find(insertRow);
  insertRowAt = insertRowAt + tmp(insertRowAt);
  channels2(insertRowAt,:) = channels(insertRow, :);
  channels2(isnan(channels2)) = channels(:);
  tIO2 = NaN(size(uniqueTime, 1) + sum(insertRow), size(uniqueTime, 2));
  tIO2(insertRowAt,:) = uniqueTime(find(insertRow)+1, :);
  tIO2(setdiff(1:(size(uniqueTime, 1) + sum(insertRow)), insertRowAt),:) = uniqueTime;

  % insert point before start of signal high (per channel)
  diffChannels = cat(1, zeros(1, size(channels2, 2), size(channels, 3)), diff(channels2));
  insertRow = any(any(diffChannels>0, 2), 3);
  channels = NaN(size(channels2, 1) + sum(insertRow), size(channels2, 2), size(channels2, 3));
  tmp = cumsum(insertRow);
  insertRowAt = find(insertRow);
  insertRowAt = insertRowAt + tmp(insertRowAt)-1;
  channels(insertRowAt,:) = channels2(find(insertRow)-1, :); % inserts "0"s
  channels(isnan(channels)) = channels2(:);
  uniqueTime = NaN(size(tIO2, 1) + sum(insertRow), size(tIO2, 2));
  uniqueTime(insertRowAt,:) = tIO2(insertRow, :);
  uniqueTime(setdiff(1:(size(tIO2, 1) + sum(insertRow)), insertRowAt),:) = tIO2;
  % make 2d matrices of correct size
  uniqueTime = repmat(uniqueTime, 1, 6);
  noChannel = all(all(channels==0, 1), 2);
  channels(:,:,noChannel) = NaN; % remove channels with no action
  channels = reshape(channels, size(channels, 1), []);

  % insert zero-level for each channel
  zeroLevel = reshape(repmat(meshgrid(1:6, 1:size(channels, 1)), size(uniqueDigiOut, 2), 1), size(channels, 1), []);
  channels(channels == 0) = zeroLevel(channels == 0)-0.5;
  channels = channels + 0.25;

  hLines = findobj(get(hAxesDigiOut, 'Children'), 'flat', 'Tag', 'plotSeqLineDigiIO');
  hPatches = findobj(get(hAxesDigiOut, 'Children'), 'flat', 'Tag', 'plotSeqPatchDigiIO');

  numLines = numel(hLines);
  numPatches = numel(hPatches);
  channelNums = 1:6;
  for iDigiOut = 1:size(uniqueTime, 2)
    % draw all unique lines
    if numLines < iDigiOut
      hLines(iDigiOut) = line(uniqueTime(:,iDigiOut), channels(:,iDigiOut), ...
        'Parent', hAxesDigiOut, 'Color', plotSequence.DigitalIOColor, ...
        'Tag', 'plotSeqLineDigiIO');
    else
      set(hLines(iDigiOut), 'XData', uniqueTime(:,iDigiOut), 'YData', channels(:,iDigiOut));
    end
  end
  iDigiOut = 0;
  for curChannel = channelNums(~noChannel)
    % draw transparent patches above channels with signals
    iDigiOut = iDigiOut +1;
    if numPatches < iDigiOut
      hPatches(iDigiOut) = patch('Parent', hAxesDigiOut, 'Faces', [1 2 3 4], ...
        'Vertices', [xLim(1), curChannel-0.25; xLim(1), curChannel+0.25; xLim(2), curChannel+0.25; xLim(2), curChannel-0.25], ...
        'FaceColor', plotSequence.DigitalIOColor, 'FaceAlpha', 0.2, ...
        'LineStyle', 'none', 'Tag', 'plotSeqPatchDigiIO', 'HitTest', 'off');
    else
      set(hPatches(iDigiOut), 'Vertices', [xLim(1), curChannel-0.25; xLim(1), curChannel+0.25; xLim(2), curChannel+0.25; xLim(2), curChannel-0.25]);
    end
  end
  % sort patches to background
  hOthers = findobj(get(hAxesDigiOut, 'Children'), 'flat', '-not', 'Tag', 'plotSeqPatchDigiIO');
  hPatches = hPatches(:);
  set(hAxesDigiOut, 'Children', [hOthers(:); hPatches(ishghandle(hPatches))]);

  % delete lines and patches that aren't used
  delete(hPatches(sum(~noChannel)+1:end));
  delete(hLines(size(uniqueTime, 2)+1:end));
end

if PlotAQ
  tempmax = 0;
  for tt = 1:length(AQ)
    t = ones(size(AQ(tt).Start,1),1) * cumsum([0,Seq.tRep(1:end-1)]);
    AQ(tt).Dur = AQ(tt).nSamples./AQ(tt).fSample;
    AQ(tt).tAQPlot = tRepStart + ...
      reshape([t(:).'+AQ(tt).Start(:).'; t(:).'+AQ(tt).Start(:).'; t(:).'+AQ(tt).Start(:).'+AQ(tt).Dur(:).'; t(:).'+AQ(tt).Start(:).'+AQ(tt).Dur(:).'], ...
      [size(AQ(tt).Start,1)*4, size(AQ(tt).Start,2), size(AQ(tt).Start,3)]);
    AQ(tt).aAQPlot = (1/HW.RX(plotSequence.iDevice).fSampleUnitScale) * ...
      reshape([NaN(size(AQ(tt).Start(:),1),1).'; AQ(tt).fSample(:).'; AQ(tt).fSample(:).'; NaN(size(AQ(tt).Start(:),1),1).'], ...
      [size(AQ(tt).Start,1)*4, size(AQ(tt).Start,2), size(AQ(tt).Start,3)]);

    tAQ = round(bsxfun(@minus, reshape(AQ(tt).tAQPlot, [], plotSequence.wraps), tStartWraps)*HW.RX(plotSequence.iDevice).fSample)/HW.RX(plotSequence.iDevice).fSample;
    fSampleAQ = reshape(AQ(tt).aAQPlot, [], plotSequence.wraps);
    uniqueTime = tAQ; uniqueTime(isnan(uniqueTime)) = 1/eps;
    uniqueFSample = fSampleAQ; uniqueFSample(isnan(uniqueFSample)) = 1/eps;

    if all(all(abs(bsxfun(@minus, uniqueTime(:,1), uniqueTime)) < 3*eps)) && ...
        all(all(abs(bsxfun(@minus, uniqueFSample(:,1), uniqueFSample)) < 3*eps))
      equalWraps = true;
      uniqueTime = tAQ(:,1);
      uniqueFSample = fSampleAQ(:,1);
      iA = 1;
    else
      equalWraps = false;

      roundTime = tAQ - min(tAQ(:));
      roundTime = uint32(roundTime/wrapDuration*(2^32-1));
      roundFSample = real(fSampleAQ) - min(fSampleAQ(:));
      roundFSample = uint32(roundFSample/max(roundFSample(:))*(2^32-1));
      [~, iA] = unique([roundTime; roundFSample].', 'rows');

      uniqueTime = tAQ(:,iA);
      uniqueFSample = fSampleAQ(:,iA);
    end

    uniqueTime = uniqueTime - plotSequence.tOffset;

    % find lines in axes
    hKids = get(hAxes(nGrads+2), 'Children');
    hLines = flipud(findobj(hKids, 'flat', 'Type', 'line', '-not', 'Tag', 'plotSeqLineLegend'));
    % compare number of found lines with required lines
    numLinesReq = size(uniqueFSample, 2);
    numLinesPre = numel(hLines)-numLinesTX*plotSequence.stackTXRX;
    if numLinesPre > numLinesReq
      % delete surplus lines
      delete(hLines(numLinesTX*plotSequence.stackTXRX+numLinesReq+1:end));
      hLines(numLinesTX*plotSequence.stackTXRX+numLinesReq+1:end) = [];
      numLinesPre = numLinesReq;
    end
    uniqueTimePre = mat2cell(uniqueTime(:,1:numLinesPre), size(uniqueTime, 1), ones(1, numLinesPre)).';
    uniqueFSamplePre = mat2cell(uniqueFSample(:,1:numLinesPre), size(uniqueFSample, 1), ones(1, numLinesPre)).';
    uniqueZeros = mat2cell(0*uniqueFSample(:,1:numLinesPre), size(uniqueFSample, 1), ones(1, numLinesPre)).';
    if plotSequence.stackTXRX
      set(hLines(numLinesTX*plotSequence.stackTXRX+1:end), {'XData'}, uniqueTimePre, {'YData'}, uniqueZeros, {'ZData'}, uniqueFSamplePre);
      % add missing lines
      hLines = [hLines; ...
        line(uniqueTime(:,numLinesPre+1:end), 0*uniqueTime(:,numLinesPre+1:end), uniqueFSample(:,numLinesPre+1:end), ...
        'Marker', '*', 'MarkerSize', 2, 'Parent', hAxes(nGrads+2))]; %#ok<AGROW>
    else
      set(hLines(numLinesTX*plotSequence.stackTXRX+1:end), {'XData'}, uniqueTimePre, {'YData'}, uniqueFSamplePre);
      % add missing lines
      hLines = [hLines; ...
        line(uniqueTime(:,numLinesPre+1:end), uniqueFSample(:,numLinesPre+1:end), ...
        'Marker', '*', 'MarkerSize', 2, 'Parent', hAxes(nGrads+2))]; %#ok<AGROW>
    end

    set(hLines(numLinesTX*plotSequence.stackTXRX+1:end), ...
      'Color', plotSequence.AQColors{tt}, 'LineWidth', 3, 'Tag', 'plotSeqLineAQ');
    % add userdata with number of samples
    tRepNumber = reshape(1:numel(Seq.tRep), [], plotSequence.wraps).';
    ws = warning('off', 'MATLAB:mat2cell:TrailingUnityVectorArgRemoved');
    [warnStr, warnId] = lastwarn();
    tStart = cellfun(@(x,y) round((x - y - plotSequence.tOffset)*HW.RX(plotSequence.iDevice).fSample)/HW.RX(plotSequence.iDevice).fSample, ...
      mat2cell(tRepStart + reshape(t(:,tRepNumber(iA,:))+AQ(tt).Start(:,tRepNumber(iA,:)), [size(AQ(tt).Start, 1), size(tRepNumber(iA,:))]), ...
      size(AQ(tt).Start, 1), ones(1, numel(iA)), size(tRepNumber, 2)), ...
      mat2cell(reshape(repmat(tStartWraps(iA), [size(AQ(tt).Start, 1), 1, size(tRepNumber(iA,:), 2)]), [size(AQ(tt).Start, 1), size(tRepNumber(iA,:))]), ...
      size(AQ(tt).Start, 1), ones(1, numel(iA)), size(tRepNumber, 2)), 'UniformOutput', false);
    nSamples = mat2cell(reshape(AQ(tt).nSamples(:,tRepNumber(iA,:)), [size(AQ(tt).nSamples, 1), size(tRepNumber(iA,:))]), ...
      size(AQ(tt).nSamples, 1), ones(1, numel(iA)), size(tRepNumber, 2));
    if isfield(AQ(tt), 'SamplingFactor')
      SamplingFactor = mat2cell(reshape(AQ(tt).SamplingFactor(:,tRepNumber(iA,:)), [size(AQ(tt).SamplingFactor, 1), size(tRepNumber(iA,:))]), ...
        size(AQ(tt).SamplingFactor, 1), ones(1, numel(iA)), size(tRepNumber, 2));
    end
    AQ_Frequency = mat2cell(reshape(AQ(tt).Frequency(:,tRepNumber(iA,:)), [size(AQ(tt).Frequency, 1), size(tRepNumber(iA,:))]), ...
      size(AQ(tt).Frequency, 1), ones(1, numel(iA)), size(tRepNumber, 2));
    AQ_Phase = mat2cell(reshape(AQ(tt).Phase(:,tRepNumber(iA,:)), [size(AQ(tt).Phase, 1), size(tRepNumber(iA,:))]), ...
      size(AQ(tt).Phase, 1), ones(1, numel(iA)), size(tRepNumber, 2));
    if ~isemptyfield(AQ(tt), 'FrequencyX') && ~isemptyfield(AQ(tt), 'PhaseX')
      AQ_FrequencyX = mat2cell(reshape(AQ(tt).FrequencyX(:,tRepNumber(iA,:)), [size(AQ(tt).FrequencyX, 1), size(tRepNumber(iA,:))]), ...
        size(AQ(tt).FrequencyX, 1), ones(1, numel(iA)), size(tRepNumber, 2));
      AQ_PhaseX = mat2cell(reshape(AQ(tt).PhaseX(:,tRepNumber(iA,:)), [size(AQ(tt).PhaseX, 1), size(tRepNumber(iA,:))]), ...
        size(AQ(tt).PhaseX, 1), ones(1, numel(iA)), size(tRepNumber, 2));
    end
    warning(ws);
    lastwarn(warnStr, warnId);
    clear ud
    ud(numLinesReq) = struct();
    [ud.tStart] = tStart{:};
    [ud.nSamples] = nSamples{:};
    if isfield(AQ(tt), 'SamplingFactor')
      [ud.SamplingFactor] = SamplingFactor{:};
    end
    [ud.Frequency] = AQ_Frequency{:};
    [ud.Phase] = AQ_Phase{:};
    if ~isemptyfield(AQ(tt), 'FrequencyX') && ~isemptyfield(AQ(tt), 'PhaseX')
      [ud.FrequencyX] = AQ_FrequencyX{:};
      [ud.PhaseX] = AQ_PhaseX{:};
    end
    ud = num2cell(ud).';
    if equalWraps
      set(hLines(numLinesTX*plotSequence.stackTXRX+1:end), 'UserData', ud{1});
    else
      set(hLines(numLinesTX*plotSequence.stackTXRX+1:end), {'UserData'}, ud);
    end

    tempmax = max(tempmax, max(AQ(tt).aAQPlot(:)));
  end
  if tempmax ~= 0 && ~plotSequence.stackTXRX
    set(hAxes(nGrads+2), 'YLim', [0, tempmax*1.1]);
  end
end

hdcm = datacursormode(hFigure);
if doRearrange && ~isempty(lastValidAxes) && ishghandle(lastValidAxes, 'axes')
  set(lastValidAxes, 'XTickLabelMode', 'auto');
  xlabel(lastValidAxes, 'time in s');
end

if ~reuseAxes
  zoom(hFigure, 'on');

  DataCursorUpdateFcns = getappdata(hFigure, 'DataCursorUpdateFcns');
  DataCursorUpdateFcns.plotSeq = @(pointDataTip, eventData) plotSeq_DataCursorFcn(pointDataTip, eventData, HW, plotSequence.iDevice);
  setappdata(hFigure, 'DataCursorUpdateFcns', DataCursorUpdateFcns);
  set(hdcm, 'UpdateFcn', @DataCursorUpdateFcnHandler);

  set(hParent, 'DeleteFcn', @plotSeq_DeleteFcn);
end
% update text in any existing data tips
hdcm.updateDataCursors();
if doRearrange && any(ishghandle(hAxes, 'axes'))
  drawnow('expose'); % necessary to have updated axes positions below (FIXME: Can this be avoided?)
  % set(hax, 'Motion', 'horizontal', 'Enable', 'on');
  % Decrease vertical spacing between axes
  firstValidAxes = hAxes(find(ishghandle(hAxes, 'axes'), 1, 'first'));
  pos = get(firstValidAxes, 'Position');
  top = pos(2) + pos(4);
  pos = get(lastValidAxes, 'Position');
  bottom = pos(2);
  spacerPart = 0.1; % portion of axes height that is used as spacer between axes
  axesHeight = (top-bottom) / (numAxes + (numAxes-1)*spacerPart);
  axesPos = get(firstValidAxes(1), 'Position');
  for iAxes = 1:numAxes
    curAxes = hAxes(axesPosition == iAxes);
    axesPos(4) = axesHeight;
    axesPos(2) = (numAxes-iAxes) * (1+spacerPart)*axesHeight + bottom;
    set(curAxes(1), 'Position', axesPos);
  end
end

% reset default zoom
if ~isequal(xLim, oldXLim) || any(newAxes)
  % FIXME: There should be a better way to set the default zoom limits...
  set(lastValidAxes, 'XLim', xLim);
  if ~isequal(xLim, oldXLim)
    hAxesResetZoom = hAxes(requiredAxes);
  else
    hAxesResetZoom = hAxes(newAxes);
  end
  for hAxis = hAxesResetZoom
    set(hFigure, 'CurrentAxes', hAxis);
    zoom(hFigure, 'reset');
  end
end

% restore previous xlims
if isequal(xLim, oldXLim) && ~isempty(prevXLim)
  set(lastValidAxes, 'XLim', prevXLim);
end

% drawnow
end


function rgbVal = str2rgb(inStr)
%% Convert character to RGB color
%
%     rgbVal = str2rgb(inStr)
%
% INPUT:
%   inStr   String or cell array of string with any of the following characters:
%           r, g, b, w, c, m, y, k.
%           Only the first character in each string (in every cell) is
%           interpreted.
%
% OUTPUT:
%   rgbVal  1x3 vector or cell array of 1x3 vectors with the corresponding RGB
%           values (range from 0 to 1).
%

colorSpec = 'rgbwcmyk';
rgbSpec = [1 0 0; 0 1 0; 0 0 1; 1 1 1; 0 1 1; 1 0 1; 1 1 0; 0 0 0];

if iscell(inStr)
  rgbVal = cellfun(@str2rgb, inStr, 'UniformOutput', false);
elseif ischar(inStr)
  colorIdx = colorSpec == inStr(1);
  if isempty(colorIdx)
    error('PD:plotSeq:WrongColor', 'Color code "%s" unknown.', inStr);
  end
  rgbVal = rgbSpec(colorIdx, :);
else
  rgbVal = inStr;
end

end


function A = fillNaNs(A, dim)
%% Replace NaNs along dimension dim with previous non-NaN value.
%
%       A = fillNaNs(A, dim)

% FIXME: Consider replacing this function by "fillmissing" once we require at
% least Matlab R2016b.

% permute dimension to operate on to front
permOrder = [dim, setdiff(1:ndims(A), dim)];
A = permute(A, permOrder);

% replace NaNs
[szDim, szOtherDims] = size(A);
for iOtherDims = 1:szOtherDims
  prevVal = A(1,iOtherDims);
  for iDim = 2:szDim
    if isnan(A(iDim,iOtherDims))
      A(iDim,iOtherDims) = prevVal;
    else
      prevVal = A(iDim,iOtherDims);
    end
  end
end

% permute dimension back to original order
A = ipermute(A, permOrder);

end


function plotSeq_DeleteFcn(hParent, eventData)
%% Executes on deletion of parent.
%
%       plotSeq_DeleteFcn(hParent, eventData)
%
% This function un-registers the DataCursorFcn when the graphics parent is
% deleted.

% FIXME: handle multiple plotSeq in one figure
hFigure = ancestor(hParent, 'figure');
DataCursorUpdateFcns = getappdata(hFigure, 'DataCursorUpdateFcns');
if isempty(DataCursorUpdateFcns)
  % FIXME: The figure ancestor might have changed since plotSeq was called. Take
  % care to carry on our own DataCursorUpdateFcn.
  return;
end
DataCursorUpdateFcns = rmfield(DataCursorUpdateFcns, 'plotSeq');
setappdata(hFigure, 'DataCursorUpdateFcns', DataCursorUpdateFcns);

end
