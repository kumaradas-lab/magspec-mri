%% Gradient Echo 2D (Flash 2D)
% This example acquires a 2D gradient echo image with Ernst angle excitation.
% 2d-image parallel to the x-z-plane of the magnet coordinate system that can be
% used to localize a sample.

%%
LoadSystem;                                         % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.Loops = 1;                                      % number of loop averages 1...

Seq.T1 = 3000e-3;                                   % T1 of sample; excitation angle is acos(exp(-Seq.tRep/Seq.T1))/pi*180
Seq.tEcho = 9e-3;                                   % echo time in seconds e.g. 4e-3
Seq.tRep = 20e-3;                                   % repetition time in seconds (default is Seq.tEcho*2)

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).Resolution = 0.1e-3;                 % resolution of image
% Seq.AQSlice(1).nRead = 16;                          % number of pixels in read direction
% Seq.AQSlice(1).nPhase(2) = 16;                      % number of pixels in phase direction
Seq.AQSlice(1).HzPerPixMin = 0;                     % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ window, 0: longest possible)
Seq.AQSlice(1).sizeRead = 0.016;                    % size in read direction in meter
Seq.AQSlice(1).sizePhase(1) = Inf;                  % size in phase(1) direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.016;                % size in phase(2) direction in meter
Seq.AQSlice(1).sizePhase(3) = Inf;                  % size in phase(3) direction in meter
Seq.AQSlice(1).thickness = 0.01;                    % slice thickness in meter
Seq.AQSlice(1).excitationPulse = @Pulse_Sinc_3_Hamming;  % excitation pulse function (type "Pulse_" than press tab for selection of pulses)

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(2) = 2;                      % oversampling phase(2)  1...

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xz');  % Top-down image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zx');  % Bottom-up image (slice/phase(1), phase(2), read/phase(3))

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeqAQ = 1:3;                                % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                   % plot every loop

Seq.AQSlice(1).plotkSpace = 1;                      % plot k-space
Seq.AQSlice(1).plotImage = 1;                       % plot image
Seq.AQSlice(1).plotImageHandle = 12;
Seq.AQSlice(1).plotPhase = 0;                       % plot phase of k-space and/or image
Seq.AQSlice(1).ZeroFillWindowSize = 2;              % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 4;                  % zero fill resolution factor
Seq.AQSlice(1).DephaseLengthFactor = 0.8;           % length factor for the dephase gradient block
Seq.AQSlice(1).SpoilLengthFactor = 1.2;             % length factor for the spoiler/crusher gradient block
Seq.AQSlice(1).sizePhaseSpoil = [1, 1, 1]*0.3e-3;   % size of spoilers/crushers in m (slice/phase(1), phase(2), read/phase(3))

Seq.CorrectSliceRephase = 0;                        % correct SliceGradTimeIntegralOffset
Seq.CorrectReadRephase = 0;                         % correct ReadGradTimeIntegralOffset
Seq.CorrectPhase = 0;                               % correct frequency drift

Seq.WeightingExponent = 1;                          % weighting exponent: 0 -> RoI weighting, 1 -> AMP.*RoI Weighting, 2 -> AMP.^2.*RoI weighting

[SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);


clf(SeqLoop.AQSlice(1).plotImageHandle);
[SeqLoop.data, SeqLoop.AQSlice] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
nR = SeqLoop.AQSlice.nRead * SeqLoop.AQSlice.ReadOS;
nP = SeqLoop.AQSlice.nPhase(2) * SeqLoop.AQSlice.PhaseOS(2);
if SeqLoop.AQSlice(1).plotkSpace
  hold(SeqLoop.AQSlice(1).plotImagehAxes{1}, 'on');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{1}, [0;nR], [floor(nP/2)+1-floor(nP/2)/2;floor(nP/2)+1-floor(nP/2)/2], '--m');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{1}, [0;nR], [floor(nP/2)+1;floor(nP/2)+1], '--m');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{1}, [0;nR], [floor(nP/2)+1+floor(nP/2)/2;floor(nP/2)+1+floor(nP/2)/2], '--m');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{1}, [floor(nR/2)+1+floor(nR/2)/2;floor(nR/2)+1+floor(nR/2)/2], [0;nP], '--m');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{1}, [floor(nR/2)+1;floor(nR/2)+1], [0;nP], '--m');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{1}, [floor(nR/2)+1-floor(nR/2)/2;floor(nR/2)+1-floor(nR/2)/2], [0;nP], '--m');
  hold(SeqLoop.AQSlice(1).plotImagehAxes{1}, 'off');
  set(SeqLoop.AQSlice(1).plotImagehAxes{1}, 'Layer', 'top');
end
if SeqLoop.AQSlice(1).plotImage
  hold(SeqLoop.AQSlice(1).plotImagehAxes{2}, 'on');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{2}, [-5;-5], [-5;5], '-m');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{2}, [0;0], [-5;5], '-m');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{2}, [5;5], [-5;5], '-m');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{2}, [-5;5], [-5;-5], '-m');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{2}, [-5;5], [0;0], '-m');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{2}, [-5;5], [5;5], '-m');
  grid(SeqLoop.AQSlice(1).plotImagehAxes{2}, 'on')
  set(SeqLoop.AQSlice(1).plotImagehAxes{2}, 'Layer', 'top', ...
    'GridColor', 'm', 'GridAlpha', 0.5, 'GridLineStyle', '--', ...
    'XMinorGrid', 'on', 'YMinorGrid', 'on', 'MinorGridColor', 'm', 'MinorGridAlpha', 0.5);

  % get center of area(s) with high intensity
  w = abs(squeeze(SeqLoop.data.ImageZ)).^SeqLoop.WeightingExponent.*SeqLoop.data.RoI;
  wRead = bsxfun(@times, w, reshape(SeqLoop.data.Ticks(1).ReadZ, [], 1));
  wPhase = bsxfun(@times, w, reshape(SeqLoop.data.Ticks(2).PhaseZ, 1, []));
  readMean = mean(wRead(:), 'omitnan') / mean(w(:), 'omitnan');
  phase2Mean = mean(wPhase(:), 'omitnan') / mean(w(:), 'omitnan');
  plot(SeqLoop.AQSlice(1).plotImagehAxes{2}, ...
    readMean/SeqLoop.AQSlice(1).LengthUnitScale, ...
    phase2Mean/SeqLoop.AQSlice(1).LengthUnitScale, 'om', ...
    'LineWidth', 2, 'MarkerSize', 8);
  hold(SeqLoop.AQSlice(1).plotImagehAxes{2}, 'off');
  LUS = SeqLoop.AQSlice(1).LengthUnitScale;
  LU = SeqLoop.AQSlice(1).LengthUnit;
  title(SeqLoop.AQSlice(1).plotImagehAxes{2}, ...
    sprintf('%c%s = %.3f %s, %c%s = %.3f %s', ...
    char(916), SeqLoop.AQSlice(1).PhaseCartesianAxis{2}, phase2Mean/LUS, LU, ...
    char(916), SeqLoop.AQSlice(1).ReadCartesianAxis{1}, readMean/LUS, LU));
end

%% -----------------------------------------------------------------------------
% (C) Copyright 2024-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
