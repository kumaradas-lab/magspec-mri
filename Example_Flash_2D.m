%% Gradient Echo 2D (Flash 2D)
% This example acquires a 2D gradient echo image with Ernst angle excitation.

%%
LoadSystem;                                         % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.Loops = 1;                                      % number of loop averages 1...

Seq.T1 = 100e-3;                                    % T1 of sample; excitation angle is acos(exp(-Seq.tRep/Seq.T1))/pi*180
Seq.tEcho = 3e-3;                                   % echo time in seconds e.g. 4e-3
Seq.tRep = 8e-3;                                   % repetition time in seconds (default is Seq.tEcho*2)

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 64;                          % number of pixels in read direction
Seq.AQSlice(1).nPhase(2) = 64;                      % number of pixels in phase direction

Seq.AQSlice(1).HzPerPixMin = 0;                     % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ window, 0: longest possible)
Seq.AQSlice(1).sizeRead = 0.01;                    % size in read direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.01;                % size in phase(2) direction in meter
Seq.AQSlice(1).thickness = 0.002;                   % slice thickness in meter
Seq.AQSlice(1).excitationPulse = @Pulse_Rect;       %Pulse_RaisedCos;  % excitation pulse function (type "Pulse_" than press tab for selection of pulses)

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(2) = 2;                      % oversampling phase(2)  1...

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orientation = 'zx';                                 % 'xy', 'yz', 'zx' for one of the cardinal planes (read-phase)
switch orientation
  case 'xy'
    Seq.AQSlice(1).alfa = 0.0*pi;                   % 1st rotation around x axis in RAD
    Seq.AQSlice(1).phi  = 0.5*pi;                   % 2nd rotation around y axis in RAD
    Seq.AQSlice(1).theta= 0.0*pi;                   % 3rd rotation around z axis in RAD
  case 'yz'
    Seq.AQSlice(1).alfa = 0.0*pi;                   % 1st rotation around x axis in RAD
    Seq.AQSlice(1).phi  = 0.0*pi;                   % 2nd rotation around y axis in RAD
    Seq.AQSlice(1).theta= 0.0*pi;                   % 3rd rotation around z axis in RAD
  case 'zx'
    Seq.AQSlice(1).alfa = 0.0*pi;                   % 1st rotation around x axis in RAD
    Seq.AQSlice(1).phi  = 0.0*pi;                   % 2nd rotation around y axis in RAD
    Seq.AQSlice(1).theta= -0.5*pi;                  % 3rd rotation around z axis in RAD
  otherwise
    if ~ischar(orientation)
      orientation = num2str(orientation);
    end
    error('Unknown orientation "%s"\n', orientation);
end
% Seq.AQSlice(1).alfa = 0.5*pi;                     % un-comment to exchange read and phase direction

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeqAQ = 1:3;                                % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                   % plot every loop
Seq.AQSlice(1).plotkSpace = 1;                      % plot k-space
Seq.AQSlice(1).plotImage = 1;                       % plot image
Seq.AQSlice(1).plotPhase = 1;                       % plot phase of k-space and/or image
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;            % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 4;                  % zero fill resolution factor

Seq.CorrectSliceRephase = 0;                        % Correct SliceGradTimeIntegralOffset
Seq.CorrectReadRephase = 0;                         % Correct ReadGradTimeIntegralOffset
Seq.CorrectPhase = 1;

[SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);

%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
