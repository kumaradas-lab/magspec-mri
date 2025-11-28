%% Measure a B1 map using two FLASH experiments (with Ernst angle and small excitation angle)

LoadSystem;


%% FLASH with Ernst angle

Seq.tEcho = 8e-3;
Seq.RepetitionTime = 20e-3;
T1 = 165e-3;  % T1 of the used sample in seconds for Ernst angle estimation
alpha_Ernst = acosd(exp(-Seq.RepetitionTime/T1));  % Ernst angle in degrees
Seq.FlipAngle = alpha_Ernst;

Seq.Loops = 1;

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 40;                          % number of pixels in read direction
Seq.AQSlice(1).nPhase(1) = 24;                      % number of pixels in phase(1) direction
Seq.AQSlice(1).nPhase(2) = 24;                      % number of pixels in phase(2) direction
Seq.AQSlice(1).HzPerPixMin = 0;                     % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ window, 0: longest possible)
Seq.AQSlice(1).sizeRead = 0.020;                    % size in read direction in meter (for CSI set to Inf)
Seq.AQSlice(1).sizePhase(1) = 0.012;                % size in phase(1) direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.012;                % size in phase(2) direction in meter
Seq.AQSlice(1).excitationPulse = @Pulse_Rect;       % excitation pulse function (type "Pulse_" than press tab for selection of pulses)

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(1) = 2;                      % oversampling phase(1)  1...
Seq.AQSlice(1).PhaseOS(2) = 2;                      % oversampling phase(2)  1...

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xyz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yzx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zxy');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeq = 1:3;                                  % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                   % plot every loop
Seq.AQSlice(1).plotkSpace = 0;                      % plot k-space
Seq.AQSlice(1).plotImage = 1201;                    % plot image
Seq.AQSlice(1).plotPhase = 0;                       % plot phase of k-space or image
Seq.AQSlice(1).plotB0ppm = 0;                       % plot B0 ppm (only 3D)
Seq.AQSlice(1).plotB0Hz = 0;                        % plot B0 Hz
Seq.AQSlice(1).ZeroFillWindowSize = 1.0;            % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 2;                  % zero fill resolution factor
Seq.AQSlice(1).RoiCutOffPercentile = 0.95;
Seq.AQSlice(1).RoiRelativeValue = 0.25;

Seq.CorrectPhase = 1;
Seq.CorrectPlotFrequency = 0;
Seq.CorrectPhaseDuration = 1e-3;
Seq.CorrectPhaseAQtOffset = 1e-3;
Seq.AQSlice(1).SpoilLengthFactor = 1.5;
Seq.AQSlice(1).DephaseLengthFactor = 0.5;

[SeqLoopErnst, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);

% display one slice in sliceomatic (if it is empty)
for iSl = 1:numel(SeqLoopErnst.AQSlice(1).plotImagehAxes)
  hsl = SeqLoopErnst.AQSlice(1).plotImagehAxes{iSl};
  if isa(hsl, 'sliceomatic') && ...
      isempty([hsl.GetAllSlicesPosX(); hsl.GetAllSlicesPosY(); hsl.GetAllSlicesPosZ(); hsl.GetAllIsoValues()])
    % Note: The axes labels don't necessarily correspond to the axes "orientation".
    hsl.AddSliceY(0);
  end
end


%% FLASH with small flip angle

% Seq.RepetitionTime = 40e-3;
Seq.FlipAngle = alpha_Ernst/4;
Seq.AQSlice(1).PhaseOS(1) = SeqLoopErnst.AQSlice(1).PhaseOS(1)*4;  % oversampling phase(1)  1...
% Seq.RepetitionTime = T1*log(cosd(Seq.FlipAngle));
Seq.AQSlice(1).plotImage = 1202;  % plot image

[SeqLoopSmall, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);

% display one slice in sliceomatic (if it is empty)
for iSl = 1:numel(SeqLoopSmall.AQSlice(1).plotImagehAxes)
  hsl = SeqLoopSmall.AQSlice(1).plotImagehAxes{iSl};
  if isa(hsl, 'sliceomatic') && ...
      isempty([hsl.GetAllSlicesPosX(); hsl.GetAllSlicesPosY(); hsl.GetAllSlicesPosZ(); hsl.GetAllIsoValues()])
    % Note: The axes labels don't necessarily correspond to the axes "orientation".
    hsl.AddSliceY(0);
  end
end


%% Calculation of B1 map
% Assumptions:
% 1. The excited signal for flip angles around the Ernst angle is approximately
%    independent of (small) flip angle deviations (i.e., small B1
%    inhomogeneities).
% 2. The excited signal is approximately linearly dependent to the flip angle
%    for small excitation angles. I.e., variations in the excited signal are
%    directly proportional to the B1 inhomogeneities of the rf coil.
% 3. The amplitude of the received signal is approximately linearly proportional
%    to B1 deviations (receive sensitivity map is proportional to the transmit
%    efficiency map for a transmit-receive-coil). In any case, the receive
%    sensitivity map is the same for the image with small excitation angle and
%    for the image at the Ernst angle.
% 4. All other deviations (sample homogeneity, B0, gradient linearity, eddy
%    currents, ...) have approximately the same influence on both acquired
%    images.
% Therefore, the ratio of the image amplitude with the small excitation angle to
% the image amplitude with the Ernst angle is directly proportional to the B1
% amplitude of the rf coil.

B1map = SeqLoopSmall.data.ImageSliceomaticZ ./ SeqLoopErnst.data.ImageSliceomaticZ;

% region of interest
RoI = SeqLoopErnst.data.RoI .* SeqLoopSmall.data.RoI;  % "and" operation
% RoI(abs(B1map) > 2) = NaN;

% cut some outer RoI voxels
extVoxelCut = 0.25;
if sum(SeqLoopErnst.AQSlice(1).nPhase(:) > 1) < 3
  % 1 read encoded and 2 phase encoded image dimensions
  RoiConvCut = SeqLoopErnst.data.ZeroFillFactor(SeqLoop.data.PermuteOrder) * extVoxelCut;
else
  % phase encoding for all three spatial dimensions (e.g., CSI)
  RoiConvCut = SeqLoopErnst.data.ZeroFillFactor(SeqLoop.data.PermuteOrder+1) * extVoxelCut;
end
if all(RoiConvCut > 0)
  convKernel = RaisedCosine(RoiConvCut*2+3);
  minNeighbors = sum(sum(sum(convKernel(:,:,ceil(size(convKernel, 3)/2 - RoiConvCut(3)):end))));
  RoI = 1 + 0 ./ (convn(~isnan(RoI), convKernel, 'same') >= minNeighbors);
end

% apply RoI to measured B1 map
B1map = B1map .* RoI;

% average ratio in RoI
meanB1au = mean(abs(B1map(:)), 'omitnan');

fprintf('average:              %.3g a.u.\n', meanB1au)
% Note: The standard deviation still includes the noise
fprintf('standard deviation:   %.3g a.u.\n', std(B1map(:), 'omitnan'))
fprintf('                      %.3g %%\n', std(B1map(:), 'omitnan')/mean(B1map(:), 'omitnan')*100)

hf = figure(1200);
B1_deviation = abs(B1map)/meanB1au*100;
hsl = sliceomatic(hf, B1_deviation, ...
  SeqLoopErnst.data.Ticks(1).ReadZ, SeqLoopErnst.data.Ticks(1).PhaseZ, SeqLoopErnst.data.Ticks(2).PhaseZ);
title(hsl.hAxes, 'B1 deviation map in %');
xlabel(hsl.hAxes, SeqLoopErnst.AQSlice(1).ReadCartesianAxis{1});
ylabel(hsl.hAxes, SeqLoopErnst.AQSlice(1).PhaseCartesianAxis{1});
zlabel(hsl.hAxes, SeqLoopErnst.AQSlice(1).PhaseCartesianAxis{2});
title(hsl.GetSliderX(), SeqLoopErnst.AQSlice(1).ReadCartesianAxis{1});
title(hsl.GetSliderY(), SeqLoopErnst.AQSlice(1).PhaseCartesianAxis{1});
title(hsl.GetSliderZ(), SeqLoopErnst.AQSlice(1).PhaseCartesianAxis{2});
% add one slice in sliceomatic (if it is empty)
if isempty([hsl.GetAllSlicesPosX(); hsl.GetAllSlicesPosY(); hsl.GetAllSlicesPosZ(); hsl.GetAllIsoValues()])
  % Note: The axes labels don't necessarily correspond to the axes "orientation".
  hsl.AddSliceY(0);
end


hf = figure(1203);
hsl = sliceomatic(hf, abs(SeqLoopErnst.data.ImageSliceomaticZ./B1map), ...
  SeqLoopErnst.data.Ticks(1).ReadZ, SeqLoopErnst.data.Ticks(1).PhaseZ, SeqLoopErnst.data.Ticks(2).PhaseZ);
title(hsl.hAxes, 'Image Amplitude (B1 corrected) in arbitrary units');
xlabel(hsl.hAxes, SeqLoopErnst.AQSlice(1).ReadCartesianAxis{1});
ylabel(hsl.hAxes, SeqLoopErnst.AQSlice(1).PhaseCartesianAxis{1});
zlabel(hsl.hAxes, SeqLoopErnst.AQSlice(1).PhaseCartesianAxis{2});
title(hsl.GetSliderX(), SeqLoopErnst.AQSlice(1).ReadCartesianAxis{1});
title(hsl.GetSliderY(), SeqLoopErnst.AQSlice(1).PhaseCartesianAxis{1});
title(hsl.GetSliderZ(), SeqLoopErnst.AQSlice(1).PhaseCartesianAxis{2});
% add one slice in sliceomatic (if it is empty)
if isempty([hsl.GetAllSlicesPosX(); hsl.GetAllSlicesPosY(); hsl.GetAllSlicesPosZ(); hsl.GetAllIsoValues()])
  % Note: The axes labels don't necessarily correspond to the axes "orientation".
  hsl.AddSliceY(0);
end


%% save measurement data and results to files
if false
  save B1_14mm.mat SeqLoopErnst SeqLoopSmall B1map B1_deviation;
  savefig(1200, 'B1_14mm_deviation_map');
  savefig(1203, 'B1_14mm_B1_Corrected_map');
end


%% -----------------------------------------------------------------------------
% (C) Copyright 2020-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
