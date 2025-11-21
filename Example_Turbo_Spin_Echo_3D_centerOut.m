%% Spin Echo 3D / Turbo Spin Echo 3D
% This sequence acquires a 3D Spin Echo image. The k-space is acquired starting
% from the center.
% It is possible to run a standard Spin Echo or a Turbo Spin Echo sequence and
% to use oversampling.

%%
LoadSystem;                                     % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.Find_Frequency_interval = 100;              % find magnet frequency every 100 seconds

Seq.Loops = 1;                                  % number of loop averages 1...

Seq.tEcho = 5e-3;                               % echo time in seconds e.g. 5e-3

Seq.average = 1;                                % number of averages
Seq.averageBreak = 100e-3;                      % break between averages

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 24;                      % number of pixels in read direction
Seq.AQSlice(1).nPhase(1) = 24;                  % number of pixels in phase(1) direction
Seq.AQSlice(1).nPhase(2) = 24;                  % number of pixels in phase(2) direction
Seq.AQSlice(1).HzPerPixMin = 500;               % bandwidth per pixel in Hz (1/HzPerPixMin= duration of AQ)
Seq.AQSlice(1).sizeRead = 0.01;                 % size in read direction in meter
Seq.AQSlice(1).sizePhase(1) = 0.01;             % size in phase(1) direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.01;             % size in phase(2) direction in meter
Seq.AQSlice(1).excitationPulse = @Pulse_Rect;   % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.AQSlice(1).inversionPulse = @Pulse_Rect;    % inversion pulse function

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(1) = 1;                  % oversampling factor in phase(1) direction
Seq.AQSlice(1).PhaseOS(2) = 2;                  % oversampling factor in phase(2) direction
Seq.AQSlice(1).oddEvenEchoes = 1;               % acquire separate images consisting of only odd or even Echoes respectively

% % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).TurboFactor = Seq.AQSlice(1).nPhase(2)*Seq.AQSlice(1).PhaseOS(2);  % number of image k-lines per excitation
Seq.AQSlice(1).TurboBreak = 100e-3;             % break between last echo and next excitation

Seq.AQSlice(1).kLineOrderType = 'centerOut';    % acquire the k-space starting from the center

Seq.AQSlice(1).SpoilLengthFactor = 0.5;         % length factor for gradient blocks adjacient to the readout
Seq.AQSlice(1).SliceRephaseLengthFactor = 2;    % length factor for slice rephase gradient
Seq.AQSlice(1).sizePhaseSpoil = [Inf, Inf, 0];  % spoil size 0 for the read spoiler means that it keeps the read-out amplitude

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xyz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yzx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zxy');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))

% % Plot        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeq = 1:3;                              % plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                               % plot every loop
Seq.AQSlice(1).plotImage = 1;                   % plot image
Seq.AQSlice(1).plotkSpace = 1;                  % plot k-space
Seq.AQSlice(1).plotPhase = 1;                   % plot phase of k-space or image
Seq.AQSlice(1).ZeroFillWindowSize = [2, 2, Inf, 2];  % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 2;              % zero fill resolution factor

[SeqLoop, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);

if 0
  %% Test zero fill and window
  SeqLoop.AQSlice(1).ZeroFillWindowSize = 1.4;    % zero fill window size (high k-space values are damped by a cos^2 law)
  SeqLoop.AQSlice(1).ZeroFillFactor = 1;          % zero fill resolution factor
  SeqLoop.AQSlice(1).plotImage = 1;               % plot image
  SeqLoop.AQSlice(1).plotPhase = 0;               % plot phase
  SeqLoop.AQSlice(1).plotkSpace = 0;              % plot k-space
  SeqLoop.AQSlice(1).tEcho = Seq.tEcho;           % echo time for B0 map
  SeqLoop.AQSlice(1).plotB0ppm = 0;               % plot B0 map in ppm
  SeqLoop.AQSlice(1).plotB0Hz = 0;                % plot B0 map in Hz
  SeqLoop.data.RoI=[];
  [SeqLoop.data] = get_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
  [SeqLoop.data] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
end

% display one slice in slicematic (if it is empty)
for iSl = 1:numel(SeqLoop.AQSlice(1).plotImagehAxes)
  hsl = SeqLoop.AQSlice(1).plotImagehAxes{iSl};
  if isa(hsl, 'sliceomatic') && ...
      isempty([hsl.GetAllSlicesPosX(); hsl.GetAllSlicesPosY(); hsl.GetAllSlicesPosZ(); hsl.GetAllIsoValues()])
    % Note: The axes labels don't necessarily correspond to the axes "orientation".
    hsl.AddSliceZ(0);
  end
end


%% -----------------------------------------------------------------------------
% (C) Copyright 2020-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
