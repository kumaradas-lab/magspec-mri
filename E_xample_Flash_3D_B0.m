%% B0 map with Gradient Echo 3D (Flash 3D)
% This sequence acquires a 3D B0 map using two gradient echo images with
% different echo times.

%%
LoadSystem;                                         % load system parameters (reset all variables HW Seq AQ TX Grad)

Seq.Loops = 0;                                      % number of loop averages 1...
Seq.LoopsBreak = 2.5;                               % Pause between two loop averages in seconds ([]= fast as possible)

Seq.T1 = 300e-3;                                    % T1 of sample; excitation angle is acos(exp(-Seq.tRep/Seq.T1))/pi*180
Seq.tEcho = 4e-3;                                   % echo time in seconds e.g. 4e-3
Seq.tRep = 2.5*Seq.tEcho+2e-3;                      % repetition time in seconds (default is Seq.tEcho*2)

% % B0 map settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.CorrectB0Read.Use = true;                       % correct read offset with B0 map
Seq.CorrectB0Read.Get = true;                       % get data for B0 correction
Seq.CorrectB0Read.tEchoIncr = 1e-3;                 % echo time increment for second measurement in s
Seq.CorrectB0Read.MinRelAmp = 0.1;
Seq.CorrectB0Read.MaxRelAmpDiff = 0.6;
% Seq.CorrectB0Read.MaxFreqOffset = 2000;
Seq.CorrectB0Read.ZeroFillWindowSize = 1;

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 16;                          % number of pixels in read direction 
Seq.AQSlice(1).nPhase(1) = 16;                      % number of pixels in phase(1) direction
Seq.AQSlice(1).nPhase(2) = 16;                      % number of pixels in phase(2) direction
Seq.AQSlice(1).HzPerPixMin = 500;                   % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ window, 0: longest possible)
Seq.AQSlice(1).sizeRead = 0.0128;                   % size in read direction in meter (for CSI set to 1e12)
Seq.AQSlice(1).sizePhase(1) = 0.0128;               % size in phase(1) direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.0128;               % size in phase(2) direction in meter
Seq.AQSlice(1).excitationPulse = @Pulse_Rect;       % excitation pulse function (type "Pulse_" than press tab for selection of pulses)

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(1) = 1;                      % oversampling phase(1)  1...
Seq.AQSlice(1).PhaseOS(2) = 2;                      % oversampling phase(2)  1... 

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).alfa = 0.0*pi;                       % 1st rotation around x axis in RAD
Seq.AQSlice(1).phi  = 0.0*pi;                       % 2nd rotation around y axis in RAD
Seq.AQSlice(1).theta= 0.0*pi;                       % 3rd rotation around z axis in RAD

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeqAQ = 1:3;                                % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                   % plot result at every loop
Seq.LoopPlotAverages = 1;                           % plot average at every loop
Seq.LoopSeqPlot = 1;
Seq.AQSlice(1).plotkSpace = 0;                      % plot k-space 
Seq.AQSlice(1).plotImage = 1;                       % plot image 
Seq.AQSlice(1).plotPhase = 0;                       % plot phase of k-space or image 
Seq.AQSlice(1).plotB0ppm = 0;                       % plot B0 ppm (only 3D)
Seq.AQSlice(1).plotB0Hz = 0;                        % plot B0 Hz  
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;            % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 2;                  % zero fill resolution factor


% HW.MagnetShim(3) = HW.MagnetShim(3) + 0.0015;

% Seq.CorrectPhase = 0;

[SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);

%% plot B0 map (before and after correction of read error)
% SeqLoop.dataB0.AQSlice.CorrectAmplitude = false;
% dataB0 = correct_read_B0(SeqLoop.dataB0, SeqLoop.AQSlice(1), SeqLoop.dataB0);

[dataB0, AQSlice] = get_kSpaceAndImageTicks(SeqLoop(1).dataB0, SeqLoop.dataB0.AQSlice);

hf = clf(figure(71));
% permuteOrder = SeqLoop(1).data(1).PermuteOrder;
permuteOrder = [2, 1, 3];
hsl = sliceomatic(hf, permute(dataB0.ImageZ, permuteOrder), ...
  dataB0.Ticks(1).ReadZ, dataB0.Ticks(1).PhaseZ, dataB0.Ticks(2).PhaseZ);
% set(hsl, 'CLim', [-1e-6, 1e-6]);
title(hsl.hAxes, 'B0 map (read corrected) in T');
xlabel(hsl.hAxes, [AQSlice.ReadCartesianAxis{1}, 'in m']);
ylabel(hsl.hAxes, [AQSlice.PhaseCartesianAxis{1}, 'in m']);
zlabel(hsl.hAxes, [AQSlice.PhaseCartesianAxis{2}, 'in m']);
title(hsl.GetSliderX(), AQSlice.ReadCartesianAxis{1});
title(hsl.GetSliderY(), AQSlice.PhaseCartesianAxis{1});
title(hsl.GetSliderZ(), AQSlice.PhaseCartesianAxis{2});

hf = clf(figure(72));
% permuteOrder = SeqLoop(1).data(1).PermuteOrder;
permuteOrder = [2, 1, 3];
dataB0.ImageZNoB0Corr(~dataB0.RoI) = NaN;
hsl = sliceomatic(hf, permute(dataB0.ImageZNoB0Corr, permuteOrder), ...
  dataB0.Ticks(1).ReadZ, dataB0.Ticks(1).PhaseZ, dataB0.Ticks(2).PhaseZ);
% set(hsl, 'CLim', [-1e-6, 1e-6]);
title(hsl.hAxes, 'B0 map (read not corrected) in T');
xlabel(hsl.hAxes, [AQSlice.ReadCartesianAxis{1}, 'in m']);
ylabel(hsl.hAxes, [AQSlice.PhaseCartesianAxis{1}, 'in m']);
zlabel(hsl.hAxes, [AQSlice.PhaseCartesianAxis{2}, 'in m']);
title(hsl.GetSliderX(), AQSlice.ReadCartesianAxis{1});
title(hsl.GetSliderY(), AQSlice.PhaseCartesianAxis{1});
title(hsl.GetSliderZ(), AQSlice.PhaseCartesianAxis{2});

return;

%%
HW.FindFrequencySweep.maxTime = 10;
Seq.dataB0 = SeqLoop.dataB0;
Seq.CorrectB0Read.Use = true;                       % correct read offset with B0 map
Seq.CorrectB0Read.Get = false;                      % get data for B0 correction
Seq.Loops = 1;                                      % number of loop averages 1...
Seq.AQSlice.CorrectAmplitude = 1;

Seq.AQSlice(1).nRead = 32;                          % number of pixels in read direction 
Seq.AQSlice(1).nPhase(1) = 32;                      % number of pixels in phase(1) direction
Seq.AQSlice(1).nPhase(2) = 32;                      % number of pixels in phase(2) direction

Seq.AQSlice(1).HzPerPixMin = 400;                  % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ window, 0: longest possible)

[SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);

%
hf = clf(figure(1110));
permuteOrder = SeqLoop(1).data(1).PermuteOrder;
% permuteOrder = [2, 1, 3];
hsl = sliceomatic(hf, permute(abs(SeqLoop(1).data(1).ImageZNoB0Corr), permuteOrder), ...
  SeqLoop(1).data(1).Ticks(1).ReadZ, SeqLoop(1).data(1).Ticks(1).PhaseZ, SeqLoop(1).data(1).Ticks(2).PhaseZ);
% set(hsl, 'CLim', [-1e-6, 1e-6]);
title(hsl.hAxes, 'Image Amplitude (read not corrected) in T');
xlabel(hsl.hAxes, [SeqLoop(1).AQSlice(1).ReadCartesianAxis{1}, 'in m']);
ylabel(hsl.hAxes, [SeqLoop(1).AQSlice(1).PhaseCartesianAxis{1}, 'in m']);
zlabel(hsl.hAxes, [SeqLoop(1).AQSlice(1).PhaseCartesianAxis{2}, 'in m']);
title(hsl.GetSliderX(), SeqLoop(1).AQSlice(1).ReadCartesianAxis{1});
title(hsl.GetSliderY(), SeqLoop(1).AQSlice(1).PhaseCartesianAxis{1});
title(hsl.GetSliderZ(), SeqLoop(1).AQSlice(1).PhaseCartesianAxis{2});

hf = clf(figure(1111));
permuteOrder = SeqLoop(1).data(1).PermuteOrder;
% permuteOrder = [2, 1, 3];
hsl = sliceomatic(hf, permute(abs(SeqLoop(1).data(1).ImageZ - SeqLoop(1).data(1).ImageZNoB0Corr), permuteOrder), ...
  SeqLoop(1).data(1).Ticks(1).ReadZ, SeqLoop(1).data(1).Ticks(1).PhaseZ, SeqLoop(1).data(1).Ticks(2).PhaseZ);
% set(hsl, 'CLim', [-1e-6, 1e-6]);
title(hsl.hAxes, 'Image Amplitude (difference) in T');
xlabel(hsl.hAxes, [SeqLoop(1).AQSlice(1).ReadCartesianAxis{1}, 'in m']);
ylabel(hsl.hAxes, [SeqLoop(1).AQSlice(1).PhaseCartesianAxis{1}, 'in m']);
zlabel(hsl.hAxes, [SeqLoop(1).AQSlice(1).PhaseCartesianAxis{2}, 'in m']);
title(hsl.GetSliderX(), SeqLoop(1).AQSlice(1).ReadCartesianAxis{1});
title(hsl.GetSliderY(), SeqLoop(1).AQSlice(1).PhaseCartesianAxis{1});
title(hsl.GetSliderZ(), SeqLoop(1).AQSlice(1).PhaseCartesianAxis{2});

%% -----------------------------------------------------------------------------
% (C) Copyright 2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
