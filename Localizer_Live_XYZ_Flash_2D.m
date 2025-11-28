%% Live 3D Localizer using 3 ortogonal 2D Gradient Echo images (Flash 2D)
%
% This example acquires 3 orthogonal 2D gradient echo images with Ernst angle
% excitation. 
% These three ortogonal 2D images can be used to localize the x, y and z
% position of a sample inside the magnet in a live interactive real-time
% measurement.

%%
LoadSystem;                                       % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.T1 = 110e-3;                                  % T1 of sample; excitation angle is acos(exp(-Seq.tRep/Seq.T1))/pi*180
Seq.tEcho = 3.5e-3;                               % echo time in seconds e.g. 4e-3
Seq.RepetitionTime = 8e-3;                        % repetition time in seconds (default is Seq.tEcho*2)
Seq.WeightingExponent = 1;                        % Weighting exponent 0 -> RoI weighting, 1 -> AMP.*RoI Weighting, 2 -> AMP.^2.*RoI weighting

Seq.T2star = Seq.T1;                              % T2star of sample in s;
Seq.T_C = 30.1;                                   % Temperatur of sample in °C;

Seq.Find_Frequency_interval = 0;
HW.FindFrequencyPause = max(Seq.T1/4, 0.5);

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).Resolution = 0.5e-3;               % image resolution in meter
% Seq.AQSlice(1).nRead = 16;                        % number of pixels in read direction
% Seq.AQSlice(1).nPhase(2) = 16;                    % number of pixels in phase direction
Seq.AQSlice(1).HzPerPixMin = 0;                   % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ window, 0: longest possible)
Seq.AQSlice(1).sizeRead = 0.025;                  % size in read direction in meter
Seq.AQSlice(1).sizePhase(1) = Inf;                % size in phase(1) direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.025;              % size in phase(2) direction in meter
Seq.AQSlice(1).sizePhase(3) = Inf;                % size in phase(3) direction in meter
Seq.AQSlice(1).thickness = 0.010;%+Inf;           % slice thickness in meter
Seq.AQSlice(1).excitationPulse = @Pulse_Sinc_3_Hamming;  % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.Slice(1).Pulse.MaxNumberOfSegments = 201;
% Seq.AQSlice(1).excitationPulse = @Pulse_RaisedCos;  % excitation pulse function (type "Pulse_" than press tab for selection of pulses)

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.Loops = 1;                                    % number of loop averages 1...
Seq.AQSlice(1).PhaseOS(2) = 1;                    % oversampling phase(2)  1...

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeqAQ = [];                               % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                 % plot every loop

Seq.AQSlice(1).plotkSpace = 0;                    % plot k-space
Seq.AQSlice(1).plotImage = 0;                     % plot image
Seq.AQSlice(1).plotImageOs = 0;                   % plot image
Seq.AQSlice(1).plotImageHandle = 11;              % Output figure number
Seq.AQSlice(1).plotPhase = 0;                     % plot phase of k-space and/or image
Seq.AQSlice(1).ZeroFillWindowSize = sqrt(2);      % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 4;                % zero fill resolution factor
Seq.AQSlice(1).DephaseLengthFactor = 0.5;         % length factor for the dephase gradient block
Seq.AQSlice(1).SpoilLengthFactor = 1;             % length factor for the spoiler/crusher gradient block
Seq.AQSlice(1).SpoilFactor = [1,1,2];             % spoiler resolution factor (slice or phase 1, phase 2, read or phase 3)
Seq.SteadyState_PreShots = 5;                     % number of excitation pulses before first acquisition

Seq.CorrectSliceRephase = 0;                      % correct SliceGradTimeIntegralOffset
Seq.CorrectReadRephase = 0;                       % correct ReadGradTimeIntegralOffset
Seq.CorrectPhase = 0;                             % correct frequency drift


%% Output Window preparation

newFigure = true;
if ishghandle(Seq.AQSlice(1).plotImageHandle)
  % set the CurrentCharacter to the Unicode replacement character which can
  % never be typed.
  set(Seq.AQSlice(1).plotImageHandle, ...
    'CurrentCharacter', native2unicode(uint32(65533), 'UTF-32'));
  newFigure = false;
end
hf = figure(Seq.AQSlice(1).plotImageHandle);

if newFigure
  % set(hf, 'WindowState', 'maximized')
  set(hf, 'Units', 'normalized', 'OuterPosition', [0.025, 0.025, 0.475, 0.95])  % [left bottom width height]
  set(hf, 'Name', [mfilename, ', press any key to stop measuring']);
end

readMean = NaN(1, 3);
phase2Mean = NaN(1, 3);

hKids = get(hf, 'Children');
hKids = hKids(ishghandle(hKids, 'axes'));
axesFound = false;
if numel(hKids) >= 4
  axTags = {'Localizer_axes1', 'Localizer_axes2', 'Localizer_axes3', '3D_Vectors'};
  axesFound = true;
  for iAx = numel(axTags):-1:1
    axc = findobj(hKids, 'Tag', axTags{iAx});
    if isempty(axc)
      axesFound = false;
      break;
    else
      hax(iAx) = axc;
    end
  end
  if isappdata(hf, 'phase2Mean')
    phase2Mean = getappdata(hf, 'phase2Mean');
  else
    axesFound = false;
  end
  if isappdata(hf, 'SeqLoop')
    SeqLoop = getappdata(hf, 'SeqLoop');
  else
    axesFound = false;
  end
  hb = findobj(hf, 'Type', 'uicontrol', 'Tag', 'RunButton');
  if isscalar(hb) && ishghandle(hb)
    set(hb, 'Value', true);
  else
    axesFound = false;
  end
end
firstRun = false;
if ~axesFound
  clf(hf);
  clear('SeqLoop', 'hax');
  hax(1) = subplot(2,2,1, 'Parent', hf);
  hax(2) = subplot(2,2,2, 'Parent', hf);
  hax(3) = subplot(2,2,3, 'Parent', hf);
  hax(4) = subplot(2,2,4, 'Parent', hf);
  firstRun = true;

  % keep measurement running while this button is pressed
  hb = uicontrol(hf, 'Style', 'togglebutton', 'Position', [10, 10, 60, 25], ...
    'String', 'Run', ...
    'Value', true, ...
    'Callback', @(h, e) restart_togglebutton_Callback(h, 'Localizer_Live_XYZ_Flash_2D'), ...
    'TooltipString', 'Use button to toggle measurement on/off.', ...
    'Tag', 'RunButton');
end
% arrayfun(@(x) hold(ax(x), 'on'), 1:4);

while 1
  for iOrient = 1:3
    % % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch abs(iOrient)
      case 1
        Seq.APT = '-yxz';
        Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, Seq.APT);  % imaging encoding directions (slice/phase(1), phase(2), read/phase(3))
        Seq.Swap = 0;
      case 2
        Seq.APT = 'zx';
        Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, Seq.APT);  % imaging encoding directions (slice/phase(1), phase(2), read/phase(3))
        Seq.Swap = 1;  % Swap read and phase when displaying the image
      case 3
        Seq.APT = 'xyz';
        Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, Seq.APT);  % imaging encoding directions (slice/phase(1), phase(2), read/phase(3))
        Seq.Swap = 0;
    end

    Seq.AQSlice(1).plotImagehAxes{2} = hax(iOrient);
    Seq.AQSlice(1).plotImageHandle = hf;
    Seq.AQSlice(1).plotImage = 0;  % plot image
    [SeqLoop(iOrient), mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);
    Seq.Find_Frequency_interval = 100;
    Seq.tRepInit = 10e-3;


    SeqLoop(iOrient).data = plot_Image(hax(iOrient), SeqLoop(iOrient).data, SeqLoop(iOrient).AQSlice(1), ...
      sprintf('Localizer_image%d', iOrient));
    hold(hax(iOrient), 'on');
    if firstRun
      plot(hax(iOrient), [5;5], [10;-10], '-m');
      plot(hax(iOrient), [10;10], [10;-10], '-m');
      plot(hax(iOrient), [-10;10], [-10;-10], '-m');
      plot(hax(iOrient), [-10;10], [-5;-5], '-m');
      plot(hax(iOrient), [-10;10], [0;0], '-m');
      plot(hax(iOrient), [-10;10], [5;5], '-m');
      plot(hax(iOrient), [-10;10], [10;10], '-m');
      plot(hax(iOrient), [-10;10], [10;10], '-m');
      plot(hax(iOrient), [-5;-5], [-10;10], '-m');
      plot(hax(iOrient), [-10;-10], [-10;10], '-m');
      plot(hax(iOrient), [0;0], [-10;10], '-m');
    else
      if exist('o', 'var') && iOrient <= numel(o) && ishghandle(o(iOrient))
        delete(o(iOrient));
      end
    end

    if SeqLoop(iOrient).Swap
      % set camera view point to swap read and phase directions
      view(hax(iOrient), 90, -90);
    end

    LUS = 1/SeqLoop(iOrient).AQSlice(1).LengthUnitScale;
    LU = SeqLoop(iOrient).AQSlice(1).LengthUnit;

    % get RoI of intens areas to cut the noise
    w = abs(squeeze(SeqLoop(iOrient).data.ImageZ)).^SeqLoop(iOrient).WeightingExponent.*SeqLoop(iOrient).data.RoI; 
    wRead = bsxfun(@times, w, reshape(SeqLoop(iOrient).data.Ticks(1).ReadZ, [], 1));
    wPhase = bsxfun(@times, w, reshape(SeqLoop(iOrient).data.Ticks(2).PhaseZ, 1, []));
    readMean(iOrient) = mean(wRead(:), 'omitnan') / mean(w(:), 'omitnan');
    phase2Mean(iOrient) = mean(wPhase(:), 'omitnan') / mean(w(:), 'omitnan');

    o(iOrient) = plot(hax(iOrient), ...
      readMean(iOrient)/SeqLoop(iOrient).AQSlice(1).LengthUnitScale, ...
      phase2Mean(iOrient)/SeqLoop(iOrient).AQSlice(1).LengthUnitScale, 'om', ...
      'LineWidth', 3, 'MarkerSize', 10); %#ok<SAGROW>
    ts = get(hax(iOrient), 'Title');
    title(hax(iOrient), ['\Delta' SeqLoop(iOrient).APT(end-1) ' = ' num2str(phase2Mean(iOrient)*LUS,'%.3f') ' ' LU]);
    hold(hax(iOrient), 'off');
    grid(hax(iOrient), 'on');
    set(hax(iOrient), 'Layer', 'top', ...
      'GridColor', 'm', 'GridAlpha', 0.5, 'GridLineStyle', '--', ...
      'XMinorGrid', 'on', 'YMinorGrid', 'on', 'MinorGridColor', 'm', 'MinorGridAlpha', 0.5);

    % Draw a coordinate system showing the displacement of the center of the weighted images.
    if iOrient==3 || ~firstRun
      plot3(hax(4), phase2Mean(2)*LUS, phase2Mean(1)*LUS, phase2Mean(3)*LUS, 'om', ...
        'LineWidth', 3, 'MarkerSize', 10);
      hold(hax(4), 'on');
      axis(hax(4), 'square');
      set(hax(4), 'DataAspectRatio', [1 1 1]);
      m = max(max(abs(phase2Mean)*1.3*LUS), 1);

      % arrows for coordinate system
      arrow3d(hax(4), [-m,0,0], [m,0,0], 10, 'cylinder', [0.1,0.2], [12,6], [0 0 1]);  % blue  [0 0 1]  [0 0 255]
      arrow3d(hax(4), [0,-m,0], [0,m,0], 10, 'cylinder', [0.1,0.2], [12,6], [1 0 0]);  % red  [1 0 0]  [255 0 0]
      arrow3d(hax(4), [0,0,-m], [0,0,m], 10, 'cylinder', [0.1,0.2], [12,6], [0 1 0]);  % green  [0 1 0]  [0 255 0]
      % arrow which shows the displacement
      arrow3d(hax(4), [phase2Mean(2),phase2Mean(1),phase2Mean(3)]*LUS, ...
        [0,0,0], 20, 'cylinder', [0.2,0.2], [12,6], [0 0 0]);  % black  [0 0 0]  [0 0 0]

      plot3(hax(4), ...
        0, 0, 0, 'ok', ...
        [0;phase2Mean(2)]*LUS, [0;phase2Mean(1)]*LUS, [0;phase2Mean(3)]*LUS, ':x');
      lw = 5;
      plot3(hax(4), [0;LUS]*phase2Mean(2), [0;0]*phase2Mean(1), [0;0]*phase2Mean(3), 'b', 'LineWidth', lw);
      plot3(hax(4), [0;LUS]*phase2Mean(2), [LUS;LUS]*phase2Mean(1), [0;0]*phase2Mean(3), ':b', 'LineWidth', lw);
      plot3(hax(4), [0;0]*phase2Mean(2), [0;LUS]*phase2Mean(1), [0;0]*phase2Mean(3), 'r', 'LineWidth', lw);
      plot3(hax(4), [LUS;LUS]*phase2Mean(2), [0;LUS]*phase2Mean(1), [0;0]*phase2Mean(3), ':r', 'LineWidth', lw);
      plot3(hax(4), [0;0]*phase2Mean(2), [0;0]*phase2Mean(1), [0;LUS]*phase2Mean(3), 'g', 'LineWidth', lw);
      plot3(hax(4), [LUS;LUS]*phase2Mean(2), [LUS;LUS]*phase2Mean(1), [0;LUS]*phase2Mean(3), ':g', 'LineWidth', lw);

      hold(hax(4), 'off');
      xlim(hax(4), [-m,m]);
      ylim(hax(4), [-m,m]);
      zlim(hax(4), [-m,m]);
      grid(hax(4), 'on');
      ylabel(hax(4), {['\Delta' SeqLoop(1).APT(end-1) ' =    '], [num2str(phase2Mean(1)*LUS,'%.3f') ' ' LU]})
      set(get(hax(4), 'YLabel'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle')
      zlabel(hax(4), {['\Delta' SeqLoop(3).APT(end-1) ' ='], [num2str(phase2Mean(3)*LUS,'%.3f') ' ' LU]})
      xlabel(hax(4), {['     \Delta' SeqLoop(2).APT(end-1) ' ='], [num2str(phase2Mean(2)*LUS,'%.3f') ' ' LU]})
      set(get(hax(4), 'XLabel'), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
      drawnow();
      
      setappdata(hf, 'phase2Mean', phase2Mean);
      setappdata(hf, 'SeqLoop', SeqLoop);
      firstRun = false;
    end

    arrayfun(@(x) set(hax(x), 'Tag', sprintf('Localizer_axes%d', x)), 1:3);
    set(hax(4), 'Tag', '3D_Vectors');
    
    % signal amplitude estimation
    if firstRun && iOrient==1
      density_H2O = get_DensityOfWaterAtTemperature(SeqLoop(iOrient).T_C);  % kg/m^3
      mol_1H = get_mol_1HAtomsFromWaterMass(density_H2O*1000);  % mol(1H)/m^3
      B1p = get_TotalMagneticPolarization(HW, mol_1H, SeqLoop(iOrient).T_C, SeqLoop(iOrient).HW.B0, SeqLoop(iOrient).HW.Gamma.H1, 0.5);

      SignalAmplitudeTheoretically = B1p ...
        * (sind(SeqLoop(iOrient).FlipAngle)*(1-exp(-(SeqLoop(iOrient).RepetitionTime)/SeqLoop(iOrient).T1))) ...
        / (1-(cosd(SeqLoop(iOrient).FlipAngle)*exp(-(SeqLoop(iOrient).RepetitionTime)/SeqLoop(iOrient).T1))) ...
        * exp(-SeqLoop(iOrient).tEcho/SeqLoop(iOrient).T2star);

      disp(['pure water signal amplitude theoretically ' num2str(SignalAmplitudeTheoretically*1e9,3) ' nT'] )
    end

    if  ((~isempty(hf.CurrentCharacter) ...
          && (hf.CurrentCharacter ~= native2unicode(uint32(65533), 'UTF-32'))) ...
         || ~get(hb, 'Value')) ...
        && ~firstRun
      break;  % break on key press for loop
    end

  end

  if  ((~isempty(hf.CurrentCharacter) ...
        && (hf.CurrentCharacter ~= native2unicode(uint32(65533), 'UTF-32'))) ...
       || ~get(hb, 'Value')) ...
      && ~firstRun
    set(hb, 'Value', false, 'Enable', 'on');
    break;  % break on key press for while loop
  end
end

%% -----------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
