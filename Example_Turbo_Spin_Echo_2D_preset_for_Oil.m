%% Spin Echo 2D / Turbo Spin Echo 2D
%
% This sequence can be used to acquire a 2D Turbo Spin Echo image.
% The echo trains are segmented in different turbo blocks to account for the
% short T2.
% The parameters are optimized for a vegetable oil sample.

%%
LoadSystem;                                                 % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.Find_Frequency_interval = 1000;                         % run Find_Frequency every X seconds or after every loop (if measurement duration > X sec)

Seq.Loops = 1;                                              % number of loop averages 1...

Seq.tEcho = 6e-3;                                           % echo time in seconds e.g. 5e-3

% % Segmentation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.T1Estimated = 200e-3;                                   % estimated T1 time of the sample in seconds; Seq.AQSlice(1).TurboBreak and Seq.LoopsBreak default to reasonable values
Seq.T2Estimated = 120e-3;                                   % estimated T2 time of the sample in seconds; Seq.AQSlice(1).TurboFactor defaults to a reasonable value

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 48;                                  % number of pixels in read; if nRead>1, nPhase(1)=1
Seq.AQSlice(1).nPhase(1) = 1;                               % number of pixels in phase(1) (3D)
Seq.AQSlice(1).nPhase(2) = 48;                              % number of pixels in phase(2) (2D)
Seq.AQSlice(1).nPhase(3) = 1;                               % number of pixels in phase(3) (CSI)
Seq.AQSlice(1).HzPerPixMin = 500;                           % bandwidth per pixel in Hz (1/HzPerPixMin= duration of AQ window)
Seq.AQSlice(1).sizeRead = 0.010;                            % size in read direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.010;                        % size in phase(2) direction in meter
Seq.AQSlice(1).thickness = 0.005;                           % slice thickness in meter
Seq.AQSlice(1).excitationPulse = @Pulse_Sinc_1;             % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.AQSlice(1).inversionPulse = @Pulse_Rect_Composite180;   % inversion pulse function

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(1) = 1;                              % oversampling phase(1)  1...
Seq.AQSlice(1).PhaseOS(2) = 2;                              % oversampling phase(2)  1...
Seq.AQSlice(1).PhaseOS(3) = 1;                              % oversampling phase(3)  1...
Seq.AQSlice(1).PhaseOS(Seq.AQSlice(1).nPhase==1) = 1;       % reset PhaseOS(x) to 1 if nPhase(x) == 1
Seq.AQSlice(1).oddEvenEchoes = 1;                           % acquire separate images consisting of only odd or even Echoes respectively. By default the average of these images is used in plot_kSpaceAndImage.

% % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The turbo factor settings are automatically set according to the segmentation
% parameters.

% % Spoiling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spoiling is probably not necessary if Seq.AQSlice(1).oddEvenEchoes is true.
% Seq.AQSlice(1).SpoilFactor = [0, 0.00, 0.00];               % spoiler resolution factor (slice or phase 1, phase 2, read or phase 3)

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xz');  % Top-down image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zx');  % Bottom-up image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeq = 1:3;                                          % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                           % plot every loop
Seq.AQSlice(1).plotkSpace = 1;                              % plot k-space
Seq.AQSlice(1).plotImage = 1;                               % plot image
Seq.AQSlice(1).plotPhase = 1;                               % plot phase of k-space or image
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;                    % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 4;                          % zero fill resolution factor

% % Other settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.MaxGradAmpSlice = 0.1;                                  % limit slice gradient strength
Seq.CorrectSliceRephase = double(Seq.AQSlice(1).thickness<=0.015);  % correct SliceGradTimeIntegralRephaseOffset
Seq.CorrectPhaseRephase = 0;                                % correct PhaseGradTimeIntegralRephaseOffset
Seq.CorrectReadRephase = 0;                                 % correct ReadGradTimeIntegralOffset

[SeqLoop, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);  % run measurement


%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
