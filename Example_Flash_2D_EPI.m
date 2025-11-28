%% (Segmented) Echo Planar Imaging EPI 2D
% This example acquires a 2D image from gradient echo trains with Ernst angle
% excitation. It is optimized for a water sample.
% The k-space is scanned mono-polar (a.k.a. uni-polar) or bi-polar in read
% direction.
% The k-lines can be acquired in multiple EPI segments.

%%
LoadSystem;                                         % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.Loops = 1;                                      % number of loop averages 1...

Seq.T1 = 3000e-3;                                   % T1 of sample for excitation with Ernst angle: acosd(exp(-Seq.tRep/Seq.T1))
Seq.tEcho = 4e-3;                                   % mean gradient echo time in seconds e.g. 4e-3
Seq.tRep = 200e-3;                                  % repetition time in seconds (default is Seq.tEcho*2)

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 65;                          % number of pixels in read direction
Seq.AQSlice(1).nPhase(2) = 16;                      % number of pixels in phase direction
Seq.AQSlice(1).HzPerPixMin = 0;                     % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ window, 0: longest possible)
Seq.AQSlice(1).sizeRead = 0.020;                    % size in read direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.020;                % size in phase(2) direction in meter
Seq.AQSlice(1).thickness = 0.003;                   % slice thickness in meter
Seq.AQSlice(1).excitationPulse = @Pulse_RaisedCos;  % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.MaxGradAmpSlice = 0.04;                         % maximum gradient amplitude for slice selection in T/m (reduce to lessen impact of eddy currents)

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(2) = 1;                      % oversampling phase(2)  1...

% % EPI settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).EPIEchoSpacing = 1e-3;               % echo spacing in EPI segments
% Seq.AQSlice(1).nImages = 2;                         % number of images (k-lines) acquired in one gradient echo train
% Seq.AQSlice(1).EPISegments = 1;                     % number of EPI segments
Seq.AQSlice(1).EPIFactor = 4;                       % EPI acceleration factor
Seq.AQSlice(1).EPIGradient = true;                  % move encoding gradient to allow them to cancel each other
% Seq.AQSlice(1).EPIGradientPhase = true;             % override setting for phase encoders
% Seq.AQSlice(1).EPIGradientRead = true;              % override setting for read encoders
Seq.AQSlice(1).EPIReadBipolar = true;               % Use bipolar read gradients in EPI segments
Seq.AQSlice(1).oddEvenEchoes = false;               % acquire separate images consisting of only odd or even Echoes respectively

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xz');  % Top-down image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zx');  % Bottom-up image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeqAQ = 1:3;                                % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                   % plot every loop
Seq.AQSlice(1).plotkSpace = 1;                      % plot k-space
Seq.AQSlice(1).plotImage = 1;                       % plot image
Seq.AQSlice(1).plotPhase = 1;                       % plot phase of k-space and/or image
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;            % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 4;                  % zero fill resolution factor

% % Optional corrections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.CorrectSliceRephase = 0;                        % correct SliceGradTimeIntegralOffset
Seq.CorrectReadRephase = 1;                         % correct ReadGradTimeIntegralOffset
Seq.CorrectReadRephasePlot = 0;                     % plot results of read rephase correction measurement
Seq.CorrectPhase = 1;                               % acquire FID directly after FID to track frequency
Seq.CorrectPlot = 0;                                % plot results of frequency (drift) measurements

[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 1);  % find magnet frequency

[SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);

if 0
  %% Test zero fill and k-space filter window
  % SeqLoop.AQSlice(1).ZeroFillWindowSize = 1.4;      % zero fill window size (high k-space values are damped by a cos^2 law)
  % SeqLoop.AQSlice(1).ZeroFillFactor = 8;            % zero fill resolution factor
  % SeqLoop.data.RoI = [];                            % RoI has to be reset in order to recalculate
  % SeqLoop.AQSlice(1).RoiCutOffPercentile = [];      % percentile that is the reference for the ROI selection
  % SeqLoop.AQSlice(1).RoiRelativeValue = 0;          % relative value of that percentile

  [SeqLoop.data] = get_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
  [SeqLoop.data, SeqLoop.AQSlice] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
end


%% -----------------------------------------------------------------------------
% (C) Copyright 2020-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
