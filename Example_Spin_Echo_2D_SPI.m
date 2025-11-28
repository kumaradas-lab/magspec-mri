%% Spin Echo 2D SPI
% This sequence acquires a 2D Spin Echo image using single point imaging (SPI).
% This can be useful if imaging a substance with short T1 and T2.

%%
LoadSystem;                                                 % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.Loops = 1;                                              % number of loop averages 1...

Seq.tEcho = 1.7e-3;                                         % echo time in seconds e.g. 5e-3
Seq.RepetitionTime = 60e-3;                                 % repetition time in seconds e.g. 100e-3 for oil or 1 for water

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 1;                                   % number of pixels in read (1: SPI)
Seq.AQSlice(1).nPhase(2) = 16;                              % number of pixels in phase(2)
Seq.AQSlice(1).nPhase(3) = 16;                              % number of pixels in phase(2)
Seq.AQSlice(1).HzPerPixMin = 1000;                          % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ)
Seq.AQSlice(1).sizeRead = Inf;                              % size in read direction in meter (Inf: SPI)
Seq.AQSlice(1).sizePhase(2) = 0.0128;                       % size in phase(2) direction in meter
Seq.AQSlice(1).sizePhase(3) = 0.0128;                       % size in phase(2) direction in meter
Seq.AQSlice(1).thickness = Inf;                             % slice thickness in meter
Seq.AQSlice(1).excitationPulse = @Pulse_RaisedCos;          % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.AQSlice(1).inversionPulse = @Pulse_Rect_Composite180;   % inversion pulse function

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The pulse sequence doesn't use spoilers (to reach shorter echo times). That
% means that the phase oversampling factor must be at least 2 to "fold" the
% inversion pulse signal out of the image.
Seq.AQSlice(1).PhaseOS(2) = 2;  % oversampling phase(2)  1...
Seq.AQSlice(1).ReadOS = 17;     % oversampling read = number of samples
Seq.AQSlice(1).ReadOSUsedForImage = Seq.AQSlice(1).ReadOS;  % number of samples that are used for image reconstruction

% % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).TurboFactor = 1;                     % number of image k-lines per excitation


% % Settings for SPI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.SingletRep = 1;                                 % group excitation and inversion into a single tRep
Seq.DephaseBefore180 = 1;                           % move phase encoding gradient to before inversion pulse to reach shorter echo times
Seq.AQSlice(1).SpoilFactor = [0, 0, 0];             % don't use spoilers to reach shorter echo times
Seq.AQSlice(1).RephaseLengthFactor = 0;             % don't rephase spin system after echo
Seq.AQSlice(1).SliceRephaseLengthFactor = 0.6;      % shorten gradient pulses to allow shorter echo times


% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orientation = 'zx';                                 % 'xy', 'yz', 'zx' for one of the cardinal planes (read-phase)
switch orientation
  case 'xy'
    Seq.AQSlice(1).alfa = 0.0*pi;                   % 1st rotation around x axis in RAD
    Seq.AQSlice(1).phi  = 0.5*pi;                   % 2nd rotation around y axis in RAD
    Seq.AQSlice(1).theta= 0.0*pi;                   % 3rd rotation around z axis in RAD
    spoilDir = 3;
  case 'yz'
    Seq.AQSlice(1).alfa = 0.0*pi;                   % 1st rotation around x axis in RAD
    Seq.AQSlice(1).phi  = 0.0*pi;                   % 2nd rotation around y axis in RAD
    Seq.AQSlice(1).theta= 0.0*pi;                   % 3rd rotation around z axis in RAD
    spoilDir = 1;
  case 'zx'
    Seq.AQSlice(1).alfa = 0.0*pi;                   % 1st rotation around x axis in RAD
    Seq.AQSlice(1).phi  = 0.0*pi;                   % 2nd rotation around y axis in RAD
    Seq.AQSlice(1).theta= -0.5*pi;                  % 3rd rotation around z axis in RAD
    spoilDir = 2;
  otherwise
    if ~ischar(orientation)
      orientation = num2str(orientation);
    end
    error('Unknown orientation "%s"\n', orientation);
end
% Seq.AQSlice(1).alfa = 0.5*pi;                     % un-comment to exchange read and phase direction

% Use continuous gradient to "spoil" the encoding
HW.Grad(1).AmpOffsetExtra(spoilDir) = HW.Grad(1).MaxAmpSlice/5;

% % Plot        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeqAQ = 1:3;                                % plot sequence where all AQs are wrapped onto each other, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                   % plot every loop
Seq.AQSlice(1).plotkSpace = 1;                      % plot k-space
Seq.AQSlice(1).plotImage = 1;                       % plot image
Seq.AQSlice(1).plotPhase = 1;                       % plot phase of k-space or image
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;            % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 4;                  % zero fill resolution factor

% Seq.CorrectSliceRephase = double(Seq.AQSlice(1).thickness<=0.015); % correct SliceGradTimeIntegralOffset
Seq.CorrectPhaseRephase = 0;                        % correct PhaseGradTimeIntegralRephaseOffset
Seq.CorrectReadRephase = 0;                         % correct ReadGradTimeIntegralOffset
Seq.MaxGradAmpSlice = 0.05;                         % limit slice gradient strength

[SeqLoop, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);

if 0
  %% Test zero fill and window
  SeqLoop.AQSlice(1).ZeroFillWindowSize = 1.4;      % zero fill window size (high k-space values are damped by a cos^2 law)
  SeqLoop.AQSlice(1).ZeroFillFactor = 8;            % zero fill resolution factor
  SeqLoop.data.RoI = [];
  [SeqLoop.data] = get_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
  [SeqLoop.data, SeqLoop.AQSlice] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
end

%% -----------------------------------------------------------------------------
% (C) Copyright 2022 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
