function spectrumData = plot_spectrum(HW, hf, ppm, spectrumData, settings)
%% Plot spectrum and drag to shift phase
%
%   spectrumData = plot_spectrum(HW, hf, ppm, spectrumData, settings)
%
% INPUT:
%   HW              HW structure
%   hf              Handle to figure with spectrum
%   ppm             x-axis
%   spectrumData    Data with complex spectrum in T
%   settings        Structure with settings for the plot. Possible fields
%                   include:
%     waitForClose    If true, blocks execution until figure is closed.
%                     (Default: false)
%     Norm2Amp        Factor for (normalized) input data when plotting.
%                     If the data is from a dual nucleus measurement, a second
%                     element can be set that applies for the X nucleus
%                     amplitudes. If is is a scalar, the value applies to both
%                     channels.
%                     (Default: 1)
%     AmpName         Name displayed on axis in plot.
%                     (Default: 'amplitude relative')
%
% OUTPUT:
%   spectrumData    If waitForClose is true, the modified complex spectrum
%                   in T is returned when the figure is closed. You can
%                   get the current complex spectrum from the open figure
%                   with handle "hf" by the following:
%                       spectrumDataCorr = getappdata(hf, 'spectrumData');
%
% Hold mouse button and move left-right for linear phase, or up-down for
% phase offset.
% Hold "Shift" for faster changes. Hold "Ctrl" for more sensitive changes.
% If the data is from a dual nucleus measurement, the left mouse button can be
% used to manipulate the data from the primary nucleus. The right mouse button
% can be used to manipulate the data from the secondary nucleus.
% Double-click to reset phase (of both nuclei).
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
if nargin < 5
  settings = struct();
end
if islogical(settings)
  wfc = settings;
  clear settings;
  settings.waitForClose = wfc;
  clear wfc;
end

if isemptyfield(settings, 'waitForClose')
  settings.waitForClose = false;
end
if isemptyfield(settings, 'Norm2Amp')
  settings.Norm2Amp = 1;
end
settings.Norm2Amp = reshape(settings.Norm2Amp, 1, []);

if isemptyfield(settings, 'AmpName')
  settings.AmpName = 'amplitude relative';
end

%%
if settings.waitForClose
  waitfor(doPlot(HW, hf, ppm, spectrumData, settings));
else
  doPlot(HW, hf, ppm, spectrumData, settings);
end

%% nested functions to keep spectrumData up-to-date in the main function
  function hf = doPlot(HW, hf, ppm, spectrumData, settings)
    %% Initialize figure with callbacks and initial plot of data
    %
    %   hf = doPlot(HW, hf, ppm, spectrumData)
    %
    % INPUT:
    %   HW              HW structure
    %   hf              Handle to figure with spectrum
    %   ppm             x-axis
    %   spectrumData    Data with complex spectrum
    %
    % OUTPUT:
    %   hf              Handle to figure with spectrum
    %

    %% initialize figure
    hf = figure(hf);
    clf(hf);
    set(hf, 'Name', 'Spectrum - Phase Shift');
    % jhf = get(hf, 'JavaFrame');
    % jj=jhf.getFigurePanelContainer.getComponent(0);
    % jj.setToolTipText('Hold mouse button and move left-right for linear phase, or up-down for phase offset')

    %% save data in figure
    setappdata(hf, 'ppm', ppm);
    setappdata(hf, 'spectrumData', spectrumData);
    setappdata(hf, 'spectrumDataOrig', spectrumData);

    %% create uipanel for controls
    isMultiFreq = size(spectrumData, 2) > 1;
    if isMultiFreq
      huic = uipanel(hf, 'Units', 'pixels', 'Position', [1, 1, 80, 60]);
      hComplex = uicontrol(huic, 'Style', 'popupmenu', ...
          'Units', 'normalized', 'Position', [0.1, 0.5, 0.8, 0.4], ...
          'String', {'all', 'real', 'abs'}, 'Value', 1, ...
          'Callback', @showComplexCore);
      setappdata(hf, 'hComplex', hComplex);
      hCore = uicontrol(huic, 'Style', 'popupmenu', ...
        'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.4], ...
        'String', {'both', '1H', 'X'}, 'Value', 1, ...
        'Callback', @showComplexCore);
      setappdata(hf, 'hCore', hCore);
    else
      huic = uipanel(hf, 'Units', 'pixels', 'Position', [1, 1, 80, 30]);
      hComplex = uicontrol(huic, 'Style', 'popupmenu', ...
          'Units', 'normalized', 'Position', [0.1, 0.05, 0.8, 0.8], ...
          'String', {'all', 'real', 'abs'}, 'Value', 1, ...
          'Callback', @showComplexCore);
      setappdata(hf, 'hComplex', hComplex);
    end

    %% plot to axes
    % axes for ticks in Hz (in background)
    haxes(4) = axes('Parent', hf, 'Position', [0.1, 0.1, 0.83, 0.85]);

    % main axes
    haxes(1) = axes('Parent', hf, 'Position', [0.1, 0.4, 0.83, 0.55]);
    % real part of spectrum
    yyaxis(haxes(1), 'left');
%     hReal = plot(haxes(1), ppm, real(spectrumData)/HW.RX.AmplitudeUnitScale);
    hReal(1) = plot(haxes(1), ppm(:,1), real(spectrumData(:,1))*settings.Norm2Amp(1), 'DisplayName', 'real');
    if isMultiFreq
      yyaxis(haxes(1), 'right');
      hReal(2:size(spectrumData, 2)) = ...
        plot(haxes(1), ppm(:,2:end), real(spectrumData(:,2:end))*settings.Norm2Amp(min(2,end)), 'DisplayName', 'real');
    end
    hold(haxes(1), 'on');
    setappdata(hf, 'hReal', hReal);
    % absolute
    yyaxis(haxes(1), 'left');
%     hAbs = plot(haxes(1), ppm, abs(spectrumData)/HW.RX.AmplitudeUnitScale, ':', 'HitTest', 'off');
    hAbs(1) = plot(haxes(1), ppm(:,1), abs(spectrumData(:,1))*settings.Norm2Amp(1), ':', 'HitTest', 'off', 'DisplayName', 'abs');
    if isMultiFreq
      yyaxis(haxes(1), 'right');
      hAbs(2:size(spectrumData, 2)) = ...
        plot(haxes(1), ppm(:,2:end), abs(spectrumData(:,2:end))*settings.Norm2Amp(min(2,end)), ':', 'HitTest', 'off', 'DisplayName', 'abs');
    end
    has_hg2 = [100, 1] * sscanf(version, '%d.', 2) >= 804;  % Before Matlab 8.4, transparency does not work
    for iLine = 1:numel(hAbs)
      lineColor = get(hAbs(iLine), 'Color');
      transp = .8;
      if has_hg2
        lineColor(4) = transp;
      else
        lineColor = (transp*lineColor + (1-transp)*[1 1 1]);  % Background is white
      end
      set(hAbs(iLine), 'Color', lineColor);
    end
    setappdata(hf, 'hAbs', hAbs);
    % imaginary part
    yyaxis(haxes(1), 'left');
%     hImag = plot(haxes(1), ppm, imag(spectrumData)/HW.RX.AmplitudeUnitScale, ':', 'HitTest', 'off');
    hImag(1) = plot(haxes(1), ppm(:,1), imag(spectrumData(:,1))*settings.Norm2Amp(1), ':', 'HitTest', 'off', 'DisplayName', 'imag');
    if isMultiFreq
      yyaxis(haxes(1), 'right');
      hImag(2:size(spectrumData, 2)) = ...
        plot(haxes(1), ppm(:,2:end), imag(spectrumData(:,2:end))*settings.Norm2Amp(min(2,end)), ':', 'HitTest', 'off', 'DisplayName', 'imag');
    end
    for iLine = 1:numel(hAbs)
      lineColor = get(hImag(iLine), 'Color');
      transp = .8;
      if has_hg2
        lineColor(4) = transp;
      else
        lineColor = (transp*lineColor + (1-transp)*[1 1 1]);  % Background is white
      end
      set(hImag(iLine), 'Color', lineColor);
    end
    setappdata(hf, 'hImag', hImag);
    % axes
    ylims = get(haxes(1), 'YLim');
%     ylabelStr = [HW.RX.AmplitudeName, ' in ' HW.RX.AmplitudeUnit];
    ylabelStr = settings.AmpName;
    ylabel(haxes(1), ylabelStr);
    if isMultiFreq
        yyaxis(haxes(1), 'left');
        ylabel(haxes(1), ylabelStr);
        yyaxis(haxes(1), 'right');
    end
    set(haxes(1), 'XTickLabel', [], 'YLim', [-ylims(2)*.5, ylims(2)]);
    legend(haxes(1), 'show');

    % axes for integral
    haxes(2) = axes('Parent', hf, 'Position', [0.1, 0.25, 0.83, 0.15]);
    acc_spec = -sign(diff(ppm(1:2)))*(bsxfun(@rdivide, cumsum(real(spectrumData), 1), sum(real(spectrumData), 1))-.5)+.5;
    hInt = plot(haxes(2), ppm, acc_spec);
    hold(haxes(2), 'on');
    acc_spec((acc_spec>=0)&(acc_spec<=1)) = NaN;
    hIntOff = plot(haxes(2), ppm, acc_spec, 'HitTest', 'off');
    set(haxes(2), 'Color', 'none', 'YLim', [-.1 1.1], 'XTickLabel', []);
    set(hInt, 'Color', [0 0 .5], 'LineWidth', 0.2);
    set(hIntOff, 'Color', 'r', 'LineWidth', 0.2);
    setappdata(hf, 'hInt', hInt);
    setappdata(hf, 'hIntOff', hIntOff);
    ylabel(haxes(2), 'spectrum integral');

    % axes for phase
    haxes(3) = axes('Parent', hf, 'Position', [0.1, 0.1, 0.83, 0.15]);
    hPhase = plot(haxes(3), ppm, angle(spectrumData));
    labelPpm = xlabel(haxes(3), 'ppm');
    set(labelPpm, 'Units', 'normalized');
    posLabelPpm = get(labelPpm, 'Position');
    posLabelPpm(1:2) = [1.025,-0.025];
    set(labelPpm, 'Position', posLabelPpm, 'HorizontalAlignment', 'left');
    ylabel(haxes(3), 'phase in rad');
    setappdata(hf, 'hPhase', hPhase);
    set(haxes(3), 'YLim', [-pi pi]);

    linkaxes(haxes(1:3), 'x');
    set(haxes(1:3), 'XGrid', 'on', 'YGrid', 'on', 'Color', 'none');
    set(haxes, 'XDir', 'reverse');
    setappdata(hf, 'haxes', haxes);

    % Add button down functions
    set(haxes(1:3), 'ButtonDownFcn', ...
        @(hObject, eventdata) axes1_ButtonDownFcn(hObject, eventdata, hf, HW));

    % settings for axes with Hz
    xlims = get(haxes(1), 'XLim');
    set(haxes(4), 'XLim', xlims/1e6*HW.fLarmor, 'YTick', [], 'XGrid', 'on', ...
        'GridLineStyle', ':', 'XAxisLocation', 'top', 'HitTest', 'off');
    transp = .9;
    gridColor = [.2 .8 .3];
    if has_hg2
      set(haxes(4), 'GridColor', gridColor, 'GridAlpha', transp);
    else
      gridColor = (transp*gridColor + (1-transp)*[1 1 1]);  % Background is white
      set(haxes(4), 'XColor', gridColor, 'Box', 'on');
    end
    labelHz = xlabel(haxes(4), 'Hz');
    set(labelHz, 'Units', 'normalized');
    posLabelHz = get(labelHz, 'Position');
    posLabelHz(1:2) = [1.025,1];
    set(labelHz, 'Position', posLabelHz, 'HorizontalAlignment', 'left');

    hZoom = zoom(hf);
    set(hZoom, 'ActionPostCallback', @(hObject, event) spectrum_ActionPostCallback(hObject, event, haxes, HW));
    hPan = pan(hf);
    set(hPan, 'ActionPostCallback', @(hObject, event) spectrum_ActionPostCallback(hObject, event, haxes, HW));
  end

  function spectrum_ActionPostCallback(hObject, eventdata, haxes, HW)
    xlims = get(haxes(1), 'XLim');
    set(haxes(4), 'XLim', xlims/1e6*HW.fLarmor);
  end

  function axes1_ButtonDownFcn(hObject, eventdata, hf, HW)
    %% Executes on Mouse button click down
    %
    %   axes1_ButtonDownFcn(hObject, eventdata, hf)
    %
    % INPUT:
    %   hObject     handle to axes1
    %   eventdata   contains info to mouse action
    %   hf          handle to figure
    %
    % Sets WindowButtonMotionFcn and WindowButtonUpFcn for the figure.
    % Resets phase on double-click.

    if strcmpi(get(hf, 'SelectionType'), 'open') % double-click
      % reset data
      spectrumData = getappdata(hf, 'spectrumDataOrig');
      setappdata(hf, 'spectrumData', spectrumData);
      updatePlot(hf, ppm, spectrumData, HW, settings);
    else
      if strcmpi(get(hf, 'SelectionType'), 'extend') % shift+click
        sensitivity = .2;
      elseif strcmpi(get(hf, 'SelectionType'), 'alt') % ctrl+click
        sensitivity = 5;
      else % normal click
        sensitivity = 1;
      end
      % FIXME: Is it ok to manipulate the signals at the two center
      % frequencies independently?
      if eventdata.Button == 3
        secondary = true;
      else
        secondary = false;
      end
      apos = get(hObject, 'CurrentPoint'); % in axes coordinates
      clickPosition = apos(1,1:2);
      set(hf, 'WindowButtonMotionFcn', @(hObject, eventdata) ...
        spectrum_WindowButtonMotionFcn(hObject, eventdata, clickPosition, sensitivity, secondary, HW));
      set(hf, 'WindowButtonUpFcn', @(hObject, eventdata) spectrum_WindowButtonUpFcn(hObject, eventdata));
    end
  end

  function spectrum_WindowButtonUpFcn(hObject, eventdata)
    %% Executes on mouse button release
    %
    %   spectrum_WindowButtonUpFcn(hObject, eventdata)
    %
    % INPUT:
    %   hObject     handle to figure
    %   eventdata   contains info to mouse action
    %
    % Clears WindowButtonMotionFcn and WindowButtonUpFcn for the figure.

    set(hObject, 'WindowButtonMotionFcn', []);
    set(hObject, 'WindowButtonUpFcn', []);
    spectrum_WindowButtonMotionFcn([]);
  end

  function spectrum_WindowButtonMotionFcn(hObject, eventdata, clickPosition, sensitivity, secondary, HW)
    %% Executes on mouse pointer movement
    %
    %   spectrum_WindowButtonMotionFcn(hObject, eventdata, clickPosition)
    %
    % INPUT:
    %   hObject         handle to figure
    %   eventdata       contains info to mouse action
    %   clickPosition   coordinates of mouse button down in axes coordinates
    %
    % Actual correction of phase on mouse movement.
    % Initialize by calling function with one empty argument.

    %% initialize
    persistent oldpos

    %% early return to initialize oldpos
    if isempty(hObject)
      oldpos = [];
      return;
    end

    %% mouse action
    apos = get(hObject, 'CurrentPoint'); % in pixels
    if isempty(oldpos)
      oldpos = apos;
      return;
    end
    pixelMoved = apos-oldpos;
    oldpos = apos;

    spectrumData = getappdata(hObject, 'spectrumData');
    % constant phase
    spectrumData(:,secondary+1) = spectrumData(:,secondary+1)*exp(-1i*pixelMoved(2)/4000/sensitivity*2*pi);
    % linear phase
    ppm = getappdata(hObject, 'ppm');
    linearPhase = pixelMoved(1)/8000/sensitivity*2*pi * (ppm(:,secondary+1) - clickPosition(1));
    spectrumData(:,secondary+1) = spectrumData(:,secondary+1).*exp(-1i*linearPhase);
    % update data
    setappdata(hObject, 'spectrumData', spectrumData);
    updatePlot(hObject, ppm, spectrumData, HW, settings);
  end
end

function updatePlot(hObject, ppm, spectrumData, HW, settings)
%% update lines in plots with spectrumData

hReal = getappdata(hObject, 'hReal');
% set(hReal, 'YData', real(spectrumData)/HW.RX.AmplitudeUnitScale);
set(hReal, {'YData'}, ...
  num2cell(bsxfun(@times, real(spectrumData), settings.Norm2Amp), 1).');
hImag = getappdata(hObject, 'hImag');
% set(hImag, 'YData', imag(spectrumData)/HW.RX.AmplitudeUnitScale);
set(hImag, {'YData'}, ...
  num2cell(bsxfun(@times, imag(spectrumData), settings.Norm2Amp), 1).');
acc_spec = -sign(diff(ppm(1:2)))*(bsxfun(@rdivide, cumsum(real(spectrumData), 1), sum(real(spectrumData), 1))-.5)+.5;
hInt = getappdata(hObject, 'hInt');
set(hInt, {'YData'}, num2cell(acc_spec, 1).');
acc_spec((acc_spec>=0)&(acc_spec<=1)) = NaN;
hIntOff = getappdata(hObject, 'hIntOff');
set(hIntOff, {'YData'}, num2cell(acc_spec, 1).');
hPhase = getappdata(hObject, 'hPhase');
set(hPhase, {'YData'}, num2cell(angle(spectrumData), 1).');

end


function showComplexCore(hObject, event)
%% callback to change visibility of displayed lines
% Select to show only real part, only absolute value, or all of real, imaginary
% and absolute value.
% Select to show only 1H, X or both spectra if applicable.

% get graphics object handles
hf = ancestor(hObject, 'figure');

hReal = getappdata(hf, 'hReal');
hImag = getappdata(hf, 'hImag');
hAbs = getappdata(hf, 'hAbs');

hComplex = getappdata(hf, 'hComplex');
if isappdata(hf, 'hCore')
  hCore = getappdata(hf, 'hCore');
end

hInt = getappdata(hf, 'hInt');
hIntOff = getappdata(hf, 'hIntOff');

hPhase = getappdata(hf, 'hPhase');

% selector with rows: real,imag,abs; and columns: 1H,X
selector = false(3, numel(hReal));

switch get(hComplex, 'Value')
  case 1  % 'all'
    selector(:,:) = true;
    set(hAbs, 'LineStyle', ':');
  case 2  % 'real'
    selector(1,:) = true;
    set(hAbs, 'LineStyle', ':');
  case 3  % 'abs'
    selector(3,:) = true;
    set(hAbs, 'LineStyle', '-');
end

if  numel(hReal) > 1
  switch get(hCore, 'Value')
    case 1  % 'both'
      set(hInt, 'Visible', 'on');
      set(hIntOff, 'Visible', 'on');
      set(hPhase, 'Visible', 'on');
    case 2  % '1H'
      set(hInt, {'Visible'}, {'on';'off'});
      set(hIntOff, {'Visible'}, {'on';'off'});
      set(hPhase, {'Visible'}, {'on';'off'});
      selector(:,2) = false;
    case 3  % 'X'
      set(hInt, {'Visible'}, {'off';'on'});
      set(hIntOff, {'Visible'}, {'off';'on'});
      set(hPhase, {'Visible'}, {'off';'on'});
      selector(:,1) = false;
  end
end

visibility = repmat({'off'}, 3, 2);
visibility(selector) = {'on'};

% change visibility
set(hReal, {'Visible'}, visibility(1,1:numel(hReal)).');
set(hImag, {'Visible'}, visibility(2,1:numel(hImag)).');
set(hAbs, {'Visible'}, visibility(3,1:numel(hAbs)).');

end
