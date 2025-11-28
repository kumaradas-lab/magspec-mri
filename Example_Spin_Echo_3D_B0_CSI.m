%% Spin Echo 3D B0 map (CSI)
% This sequence acquires a B0 map using a CSI 3D Spin Echo measurement.

%%
LoadSystem;                                     % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.Loops = 1;                                  % number of loop averages 1...

Seq.tEcho = 5e-3;                               % Echo time in seconds e.g. 5e-3
Seq.RepetitionTime = 80e-3;                     % repetition time in seconds


% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 1;                       % number of pixels in read direction
Seq.AQSlice(1).nPhase(1) = 10;                  % number of pixels in phase(1) direction
Seq.AQSlice(1).nPhase(2) = 20;                  % number of pixels in phase(2) direction
Seq.AQSlice(1).nPhase(3) = 10;                  % number of pixels in phase(3) direction

Seq.AQSlice(1).HzPerPixMin = 500;               % bandwidth per pixel in Hz (1/HzPerPixMin= duration of AQ)
Seq.AQSlice(1).sizeRead = Inf;                  % size in read direction in meter
Seq.AQSlice(1).sizePhase(1) = 0.010;            % size in phase(1) direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.020;            % size in phase(2) direction in meter
Seq.AQSlice(1).sizePhase(3) = 0.010;            % size in phase(3) direction in meter
Seq.AQSlice(1).excitationPulse = @Pulse_Rect;   % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.AQSlice(1).inversionPulse = @Pulse_Rect;    % inversion pulse function

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(1) = 1;                  % oversampling phase(1)
Seq.AQSlice(1).PhaseOS(2) = 1;                  % oversampling phase(2)
Seq.AQSlice(1).PhaseOS(3) = 1;                  % oversampling phase(3)
Seq.AQSlice(1).ReadOS = 51;                     % oversampling read

% % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).TurboFactor = 1;                 % number of image k-lines per excitation
% Seq.AQSlice(1).TurboBreak = Seq.tRep;           % break between last echo in echo train and next excitation

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xyz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yzx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zxy');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))

Seq.SteadyState_PreShots90 = 2;
Seq.SteadyState_PostShots90 = 0;

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeq = 1:3;                              % plot pulse program, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.plotSequence.wraps = prod(Seq.AQSlice(1).nPhase(2:3))*prod(Seq.AQSlice(1).PhaseOS(2:3));
Seq.plotSeqStart = 2*Seq.SteadyState_PreShots90+1;
Seq.plotSeqEnd = 2*(prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS)+Seq.SteadyState_PreShots90);

Seq.LoopPlot = 1;                               % plot every loop
Seq.AQSlice(1).plotImage = 1;                   % plot image
Seq.AQSlice(1).plotkSpace = 0;                  % plot k-space
Seq.AQSlice(1).plotPhase = 0;                   % plot phase of k-space or image
Seq.AQSlice(1).plotB0PpmPhase = 1;              % plot B0 ppm (only 3D)
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;        % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 4;              % zero fill resolution factor
Seq.AQSlice(1).ReadOSUsedForImage = 1;          % number of (oversampled) samples used for image


[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);  % find magnet frequency

HW.Grad(1).AmpOffset(1:3) = HW.Grad(1).AmpOffset(1:3) + [0, 0, 0]*1e-3;

[SeqLoop, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);

LimitsPPM = [-10, 10]; % in ppm
set(SeqLoop.AQSlice(1).plotB0ppmhAxes{2}, 'CLim', LimitsPPM);

B0PpmPhase = SeqLoop.data.RoI.*unwrap3Dmiddle(-angle(SeqLoop.data.ImageSliceomaticZ))/pi/2/SeqLoop.tEcho/SeqLoop.data.fCenter*1e6;

% display one slice in slicematic (if it is empty)
for iSl = 1:numel(SeqLoop.AQSlice(1).plotImagehAxes)
  hsl = SeqLoop.AQSlice(1).plotImagehAxes{iSl};
  if isa(hsl, 'sliceomatic') && ...
      isempty([hsl.GetAllSlicesPosX(); hsl.GetAllSlicesPosY(); hsl.GetAllSlicesPosZ(); hsl.GetAllIsoValues()])
    % Note: The axes labels don't necessarily correspond to the axes "orientation".
    hsl.AddSliceZ(0);
  end
end
for iSl = 1:numel(SeqLoop.AQSlice(1).plotB0ppmhAxes)
  hsl = SeqLoop.AQSlice(1).plotB0ppmhAxes{iSl};
  if isa(hsl, 'sliceomatic') && ...
      isempty([hsl.GetAllSlicesPosX(); hsl.GetAllSlicesPosY(); hsl.GetAllSlicesPosZ(); hsl.GetAllIsoValues()])
    % Note: The axes labels don't necessarily correspond to the axes "orientation".
    hsl.AddSliceZ(0);
  end
end


%% -----------------------------------------------------------------------------
% (C) Copyright 2018-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
