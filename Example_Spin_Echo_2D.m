%% Spin Echo 2D / Turbo Spin Echo 2D
% This sequence acquires a 2D Spin Echo image.
% It is possible to run a standard SE or a Turbo SE and use oversampling

%%
LoadSystem;                                                 % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.Loops = 1;                                              % number of loop averages 1...

Seq.tEcho = 5e-3;                                           % echo time in seconds e.g. 5e-3
Seq.RepetitionTime = 100e-3;                                % repetition time in seconds e.g. 100e-3 for oil or 1 for water

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 64;                                  % number of pixels in read, if nRead>1 nPhase(1)=1
Seq.AQSlice(1).nPhase(2) = 64;                              % number of pixels in phase(2)
Seq.AQSlice(1).HzPerPixMin = 500;                           % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ)
Seq.AQSlice(1).sizeRead = 0.0128;                           % size in read direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.0128;                       % size in phase(2) direction in meter
Seq.AQSlice(1).thickness = 0.005;                           % slice thickness in meter
Seq.AQSlice(1).excitationPulse = @Pulse_RaisedCos;          % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.AQSlice(1).inversionPulse = @Pulse_Rect_Composite180;   % inversion pulse function

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(2) = 2;                              % oversampling phase(2)  1...

% % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).TurboFactor = 1;                             % number of image k-lines per excitation
% Seq.AQSlice(1).TurboBreak = Seq.RepetitionTime;             % break between last echo and next excitation

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xz');  % Top-down image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zx');  % Bottom-up image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))

% % Plot        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeqAQ = 1:3;                                        % plot sequence where all AQs are wrapped onto each other, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                           % plot every loop
Seq.AQSlice(1).plotkSpace = 1;                              % plot k-space
Seq.AQSlice(1).plotImage = 1;                               % plot image
Seq.AQSlice(1).plotPhase = 1;                               % plot phase of k-space or image
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;                    % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 4;                          % zero fill resolution factor

% Seq.CorrectSliceRephase = double(Seq.AQSlice(1).thickness<=0.015); % correct SliceGradTimeIntegralOffset
Seq.CorrectPhaseRephase = 0;                                % correct PhaseGradTimeIntegralRephaseOffset
Seq.CorrectReadRephase = 0;                                 % correct ReadGradTimeIntegralOffset
Seq.MaxGradAmpSlice = 0.05;                                 % limit slice gradient strength

[SeqLoop, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);

if 0
  %% Test zero fill and window
  SeqLoop.AQSlice(1).ZeroFillWindowSize = 1.4;                % zero fill window size (high k-space values are damped by a cos^2 law)
  SeqLoop.AQSlice(1).ZeroFillFactor = 8;                      % zero fill resolution factor
  SeqLoop.data.RoI = [];
  [SeqLoop.data] = get_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
  [SeqLoop.data, SeqLoop.AQSlice] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
end


%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
