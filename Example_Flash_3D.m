%% Gradient Echo 3D (Flash 3D)
% This sequence acquires a 3D gradient echo image.

%%
LoadSystem;                                     % load system parameters (reset all variables HW Seq AQ TX Grad)

Seq.Loops = 1;                                  % number of loop averages 1...

Seq.T1 = 100e-3;                                % T1 of sample; excitation angle is acos(exp(-Seq.tRep/Seq.T1))/pi*180
Seq.tEcho = 3e-3;                               % echo time in seconds e.g. 4e-3
Seq.tRep = 8e-3;                               % repetition time in seconds (default is Seq.tEcho*2)

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 32;                      % number of pixels in read direction
Seq.AQSlice(1).nPhase(1) = 32;                  % number of pixels in phase(1) direction
Seq.AQSlice(1).nPhase(2) = 32;                  % number of pixels in phase(2) direction
Seq.AQSlice(1).HzPerPixMin = 0;                 % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ window, 0: longest possible)
Seq.AQSlice(1).sizeRead = 0.01;                 % size in read direction in meter (for CSI set to 1e12)
Seq.AQSlice(1).sizePhase(1) = 0.01;             % size in phase(1) direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.01;             % size in phase(2) direction in meter
Seq.AQSlice(1).excitationPulse = @Pulse_Rect;   % excitation pulse function (type "Pulse_" than press tab for selection of pulses)

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(1) = 2;                  % oversampling phase(1)  1...
Seq.AQSlice(1).PhaseOS(2) = 2;                  % oversampling phase(2)  1...

% % Spoiler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).SpoilLengthFactor = 2;

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xyz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yzx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zxy');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeq = 1:3;                              % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                               % plot every loop
Seq.AQSlice(1).plotkSpace = 0;                  % plot k-space
Seq.AQSlice(1).plotImage = 1;                   % plot image
Seq.AQSlice(1).plotPhase = 0;                   % plot phase of k-space or image
Seq.AQSlice(1).plotB0ppm = 0;                   % plot B0 ppm (only 3D)
Seq.AQSlice(1).plotB0Hz = 0;                    % plot B0 Hz
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;        % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 2;              % zero fill resolution factor

[SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);

% display one slice in sliceomatic (if it is empty)
for iSl = 1:numel(SeqLoop.AQSlice(1).plotImagehAxes)
  hsl = SeqLoop.AQSlice(1).plotImagehAxes{iSl};
  if isa(hsl, 'sliceomatic') && ...
      isempty([hsl.GetAllSlicesPosX(); hsl.GetAllSlicesPosY(); hsl.GetAllSlicesPosZ(); hsl.GetAllIsoValues()])
    % Note: The axes labels don't necessarily correspond to the axes "orientation".
    hsl.AddSliceZ(0);
  end
end


%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
