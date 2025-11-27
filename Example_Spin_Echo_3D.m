%% Spin Echo 3D / Turbo Spin Echo 3D
% This sequence acquires a 3D Spin Echo image.
% It is possible to run a standard Spin Echo or a Turbo Spin Echo sequence and
% to use oversampling.

%%
LoadSystem;                                     % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.Loops = 1;                                  % number of loop averages 1...

Seq.tEcho = 5e-3;                               % echo time in seconds e.g. 5e-3
Seq.RepetitionTime = 100e-3;                    % repetition time in seconds

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 16;                      % number of Pixels in read direction
Seq.AQSlice(1).nPhase(1) = 16;                  % number of Pixels in phase(1) direction
Seq.AQSlice(1).nPhase(2) = 16;                  % number of Pixels in phase(2) direction
Seq.AQSlice(1).HzPerPixMin = 500;               % bandwidth per pixel in Hz (1/HzPerPixMin= duration of AQ)
Seq.AQSlice(1).sizeRead = 0.01;                 % size in read direction in meter
Seq.AQSlice(1).sizePhase(1) = 0.01;             % size in phase(1) direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.01;             % size in phase(2) direction in meter
Seq.AQSlice(1).excitationPulse = @Pulse_Rect;   % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.AQSlice(1).inversionPulse = @Pulse_Rect;    % inversion pulse function

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(1) = 1;                  % Oversampling phase(1)  1...
Seq.AQSlice(1).PhaseOS(2) = 2;                  % Oversampling phase(2)  1...

% % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).TurboFactor = 1;                 % number of image k-lines per excitation
% Seq.AQSlice(1).TurboBreak = Seq.RepetitionTime;  % break between last echo and next excitation

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).alfa = 0.0*pi;                   % 1st rotation around x axis in RAD
Seq.AQSlice(1).phi = 0.0*pi;                    % 2nd rotation around y axis in RAD
Seq.AQSlice(1).theta = 0.0*pi;                  % 3rd rotation around z axis in RAD

% % Plot        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeq = 1:3;                              % plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                               % plot every loop
Seq.AQSlice(1).plotImage = 1;                   % plot image
Seq.AQSlice(1).plotkSpace = 0;                  % plot k-space
Seq.AQSlice(1).plotPhase = 0;                   % plot phase of k-space or image
Seq.AQSlice(1).ZeroFillWindowSize = 2;          % zero fill window size (high k-space values are damped by a cos^2 law)
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
% (C) Copyright 2011-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
