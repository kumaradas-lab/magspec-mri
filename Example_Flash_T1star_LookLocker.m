%% Gradient Echo (Flash) image with inversion pulses for T1 measurements
% This example acquires multiple 2-d gradient echo images with narrow angle
% excitation after a preparing inversion pulse.
% This is essentially the method proposed by Look and Locker in 1970.
% The "apparent T1" (T1*) is determined with an exponential fit to the
% amplitudes at each pixel in the resulting images.


%%
LoadSystem;                                         % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.T1Estimated = 200e-3;                           % T1 of sample; excitation angle is acos(exp(-Seq.tRep/Seq.T1Estimated))/pi*180

Seq.Loops = 1;                                      % number of loop averages 1...

% % Look-Locker pulse sequence settings  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.Function_Prepare_Measurement = @prepare_LookLocker;
Seq.tRelax = Seq.T1Estimated*5;  % relaxation time before preparation pulses in seconds
Seq.FlipAngle = 4;  % flip angle in degrees
Seq.AQSlice(1).nImagesSlice = 5;  % number of small angle excitations after each preparation pulse
Seq.tRep = Seq.T1Estimated*3/Seq.AQSlice(1).nImagesSlice;  % repetition time of small excitation pulses in seconds
Seq.tPrepare = 10e-3;  % time between preparation pulse and first excitation pulse in seconds
Seq.tEcho = 3.5e-3;  % echo time in seconds e.g. 4e-3
Seq.preparationSpoilSize = 0.02e-3;  % size of spoiler after preparation pulse in meter
Seq.preparationSpoilLength = min(5e-3, Seq.tPrepare-2e-3);  % length of spoiler after preparation pulse in seconds
Seq.SteadyState_PreShots = 2 * Seq.AQSlice(1).nImagesSlice;  % number of preshots before actual measurement

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 24;                          % number of pixels in read direction
Seq.AQSlice(1).nPhase(1) = 1;                       % number of pixels in phase(1) direction (enables 3-d)
Seq.AQSlice(1).nPhase(2) = 24;                      % number of pixels in phase(2) direction (2-d and 3-d)
Seq.AQSlice(1).HzPerPixMin = 0;                     % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ window, 0: longest possible)
Seq.AQSlice(1).sizeRead = 0.010;                    % size in read direction in meter
Seq.AQSlice(1).sizePhase(1) = 0.010;                % size in phase(1) direction in meter (only 3-d)
Seq.AQSlice(1).sizePhase(2) = 0.010;                % size in phase(2) direction in meter (2-d and 3-d)
if sum(Seq.AQSlice(1).nPhase > 1) < 2
  % 2-d image
  Seq.AQSlice(1).thickness = 5e-3;                  % slice thickness in meter
end
Seq.AQSlice(1).excitationPulse = @Pulse_RaisedCos;  % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.AQSlice(1).SpoilLengthFactor = 10;

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(1) = 1;                      % oversampling phase(1)  1...
Seq.AQSlice(1).PhaseOS(2) = 2;                      % oversampling phase(2)  1...

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(Seq.AQSlice(1).nPhase > 1) < 2
  % 2-d image
  Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xz');  % Top-down image (phase(2), read/phase(3))
  % Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zx');  % Bottom-up image (phase(2), read/phase(3))
  % Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yx');  % image encoding directions (phase(2), read/phase(3))
  % Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yz');  % image encoding directions (phase(2), read/phase(3))
else
  % 3-d image
  Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xyz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
  % Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yzx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
  % Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zxy');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
end

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeqAQ = 1:3;                                % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                   % plot every loop
Seq.AQSlice(1).plotkSpace = 1;                      % plot k-space
Seq.AQSlice(1).plotImage = 1;                       % plot image
Seq.AQSlice(1).plotPhase = 1;                       % plot phase of k-space and/or image
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;            % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 2;                  % zero fill resolution factor

% % some corrections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.CorrectSliceRephase = 0;  % correct SliceGradTimeIntegralOffset
Seq.CorrectReadRephase = 0;  % correct ReadGradTimeIntegralOffset
Seq.CorrectPhase = 1;  % acquire frequency tracking windows to correct the phase of the acquired signal
Seq.CorrectPhaseDuration = 1e-3;  % duration of frequency tracking windows in s
Seq.CorrectPhase_maxPhaseOffset = 2*pi/25;  % maximum tolerated phase offset in rad due to standard error of frequency tracking

% % post-processing settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.PhaseCorrectionWithLastImage = true;  % use phase of last image as reference

[SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);

if 0
  %% Test zero fill and k-space filter window
  % SeqLoop.AQSlice(1).ZeroFillWindowSize = 1.4;      % zero fill window size (high k-space values are damped by a cos^2 law)
  % SeqLoop.AQSlice(1).ZeroFillFactor = 8;            % zero fill resolution factor
  % SeqLoop.data.RoI = [];                            % RoI has to be reset in order to recalculate
  % SeqLoop.AQSlice(1).RoiCutOffPercentile = [];      % percentile that is the reference for the ROI selection
  % SeqLoop.AQSlice(1).RoiRelativeValue = 0;          % relative value of that percentile
  SeqLoop.AQSlice(1).iSlice = 2;                    % number of the image in the series to show

  [SeqLoop.data] = get_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
  [SeqLoop.data, SeqLoop.AQSlice] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
end

%% post-processing
% TODO: Can probably be very similar to the post-processing in
% Example_Spin_Echo_T2.m.

clear data
data.ImageZ = SeqLoop.data.ImageZ; % get images
data.kSpaceOsRaw = SeqLoop.data.kSpaceOsRaw; % get k-space of all images

if SeqLoop.PhaseCorrectionWithLastImage
  % assume that the signal is sufficiently recovered for the last acquired image
  % in the series
  data.ImageZFirstPhaseCor = exp(-1i*(0+angle(data.ImageZ(:,:,:,:,end))));
  data.kSpaceOsRawFirstPhaseCor = exp(-1i*(0+angle(data.kSpaceOsRaw(:,:,:,:,end))));
  % phase correction with first image and k-space
  data.ImageZ = bsxfun(@times, data.ImageZ, exp(-1i*(0+angle(data.ImageZ(:,:,:,:,end)))));
  imax = find(data.kSpaceOsRaw(:,:,:,:,end) == max(max(max(max(data.kSpaceOsRaw(:,:,:,:,end))))), 1, 'first');
  [i1, i2, i3, i4] = ind2sub(size(data.kSpaceOsRaw(:,:,:,:,end)), imax);
  data.kSpaceOsRaw = data.kSpaceOsRaw * exp(-1i*(0+angle(data.kSpaceOsRaw(i1,i2,i3,i4,end))));
else
  data.ImageZFirstPhaseCor = ones(size(data.ImageZFirstPhaseCor));
  data.kSpaceOsRawFirstPhaseCor = ones(size(data.kSpaceOsRawFirstPhaseCor));
end

% permute dimension with tEchoes to front in amplitude
amax = max(abs(data.ImageZ(:)));
data.DataAmplitude = permute(data.ImageZ(:,:,:,:,1:end), [5,1,2,3,4])/amax;
data.DataTime = mean(SeqLoop.data.tImageZ(1:end,:), 2) + SeqLoop.tPrepare;

% Plot some images of the T1 encoded image stack
minmax = [-max(abs(data.ImageZ(:))), max(abs(data.ImageZ(:)))];
amaxk = max(abs(data.kSpaceOsRaw(:)));
minmaxk = [-max(abs(data.kSpaceOsRaw(:))), max(abs(data.kSpaceOsRaw(:)))];
if SeqLoop.AQSlice(1).nImagesSlice < 8
  numAll = 1:SeqLoop.AQSlice(1).nImagesSlice;
elseif SeqLoop.AQSlice(1).nImagesSlice < 16
  numAll = [1,2,3,6:3:SeqLoop.AQSlice(1).nImagesSlice];
elseif SeqLoop.AQSlice(1).nImagesSlice < 24
  numAll = [1,2,4,7:4:SeqLoop.AQSlice(1).nImagesSlice];
elseif SeqLoop.AQSlice(1).nImagesSlice < 8*5
  numAll = [1,2,4,9:6:SeqLoop.AQSlice(1).nImagesSlice];
elseif SeqLoop.AQSlice(1).nImagesSlice < 8*6
  numAll=[1,2,5,12:7:SeqLoop.AQSlice(1).nImagesSlice];
else % if SeqLoop.AQSlice(1).nImagesSlice < 57
  numAll = [1,2,6,13,20:9:SeqLoop.AQSlice(1).nImagesSlice];
end

iAxes = 0;
dim_image = sum([SeqLoop.AQSlice(1).nPhase, SeqLoop.AQSlice(1).nRead] > 1);
if dim_image < 3
  % 2d image
  if Seq.AQSlice(1).plotkSpace, k=2; else k=1; end
  hf = figure(2);
  clf(hf, 'reset')
  for num = numAll
    iAxes = iAxes+1;
    minmax2 = [-max(max(max(max(max(abs(data.ImageZ(:,:,:,:,num)))))))-eps(amax), ...
      max(max(max(max(max(abs(data.ImageZ(:,:,:,:,num)))))))+eps(amax)]/amax;
    % image amplitude
    hax = subplot(numel(numAll),5*k,1+(iAxes-1)*5*k);
    imagesc(squeeze(abs(data.ImageZ(:,:,:,:,num))).'/amax, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmax2);
    title(hax, sprintf('abs %.1f', 1/minmax2(2)));  % scale factor for amp
    ylabel(hax, sprintf('TRel%d = %.1f ms', num, SeqLoop.data.tImageZ(num)*1e3));  % relaxation time after preparation pulse
    % image real
    hax = subplot(numel(numAll),5*k,2+(iAxes-1)*5*k, 'Parent', hf);
    imagesc(squeeze(real(data.ImageZ(:,:,:,:,num))).'/amax, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmax2);
    title(hax, sprintf('real %.1f', 1/minmax2(2)));  % scale factor for amp
    % image imag
    hax = subplot(numel(numAll),5*k,3+(iAxes-1)*5*k, 'Parent', hf);
    imagesc(squeeze(imag(data.ImageZ(:,:,:,:,num))).'/amax, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmax2);
    title(hax, sprintf('imag %.1f', 1/minmax2(2)));
    % image phase
    hax = subplot(numel(numAll),5*k,4+(iAxes-1)*5*k, 'Parent', hf);
    imagesc(squeeze(angle(data.ImageZ(:,:,:,:,num))).', 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', [-pi,pi]);
    title(hax, 'phase');
    % image amplitude with color limits of first row
    hax = subplot(numel(numAll),5*k,5+(iAxes-1)*5*k, 'Parent', hf);
    imagesc(squeeze(abs(data.ImageZ(:,:,:,:,num))).'/amax, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmax/amax);
    title(hax, sprintf('abs %.1f, const. CLim', 1/minmax2(2)));
    if k==2
      minmaxk2 = [-max(max(max(max(max(abs(data.kSpaceOsRaw(:,:,:,:,num)))))))-eps(amaxk), ...
        max(max(max(max(max(abs(data.kSpaceOsRaw(:,:,:,:,num)))))))+eps(amaxk)]/amaxk;
      % kSpace amplitude
      hax = subplot(numel(numAll),5*k,1+(iAxes-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(abs(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', minmaxk2);
      title(hax, sprintf('abs %.1f', 1/minmaxk2(2)));
      % kSpace real
      hax = subplot(numel(numAll),5*k,2+(iAxes-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(real(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', minmaxk2);
      title(hax, sprintf('real %.1f', 1/minmaxk2(2)));

      hax = subplot(numel(numAll),5*k,3+(iAxes-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(imag(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', minmaxk2);
      title(hax, sprintf('imag %.1f', 1/minmaxk2(2)));

      hax = subplot(numel(numAll),5*k,4+(iAxes-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(angle(data.kSpaceOsRaw(:,:,:,:,num))).', 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', [-pi,pi]);
      title(hax, 'phase');

      hax = subplot(numel(numAll),5*k,5+(iAxes-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(abs(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', minmaxk/amaxk);
      title(hax, 'abs 1.00');
    end

  end
  drawnow();
end

% get region of interest (RoI) from last image in series (assuming magnetization
% is relaxed at that point)
if isemptyfield(SeqLoop.AQSlice(1), 'RoiRelativeValue')
  SeqLoop.AQSlice(1).RoiRelativeValue = 1/5;
end
if isemptyfield(SeqLoop.AQSlice(1), 'RoiCutOffPercentile')
  SeqLoop.AQSlice(1).RoiCutOffPercentile = 0.95;
end
if isemptyfield(SeqLoop.AQSlice(1), 'RoiFitAmpToOffsetFactorTol')
  SeqLoop.AQSlice(1).RoiFitAmpToOffsetFactorTol = 1.4;
end
imageSorted = sort(abs(reshape(SeqLoop.data.ImageZ(:,:,:,:,end), [], 1)));
SeqLoop.data.RoICutOff = imageSorted(round(numel(SeqLoop.data.ImageZ(:,:,:,:,end))*SeqLoop.AQSlice(1).RoiCutOffPercentile));
clear imageSorted
SeqLoop.data.RoIZ = ones(size(SeqLoop.data.ImageZ));
SeqLoop.data.RoIZ(repmat(abs(SeqLoop.data.ImageZ(:,:,:,:,end))<=SeqLoop.data.RoICutOff.*SeqLoop.AQSlice(1).RoiRelativeValue, ...
                         1, 1, 1, 1, size(SeqLoop.data.ImageZ, 5))) = 0;
kernel_sz = SeqLoop.data.ZeroFillFactor*2+1;
kernel_sz(SeqLoop.data.ZeroFillFactor==1) = 1;
SeqLoop.data.RoIZ = convn(SeqLoop.data.RoIZ, ones(kernel_sz), 'same');
SeqLoop.data.RoIZ(SeqLoop.data.RoIZ<prod(kernel_sz)*2/3) = NaN;
SeqLoop.data.RoIZ(~isnan(SeqLoop.data.RoIZ)) = 1;
SeqLoop.data.RoIZ = squeeze(SeqLoop.data.RoIZ);


%% Single-exponential fit of amplitude at each pixel
nImages = size(data.DataAmplitude, 1);
if dim_image < 3
  % 2d image
  RoI = SeqLoop.data.RoIZ(:,:,1);
  SeqLoop.data.RoIAll = bsxfun(@times, reshape(RoI, 1, size(RoI,1), 1, size(RoI,2)), ...
    [ones(nImages, 1); NaN(size(data.DataAmplitude,1)-nImages, 1)]);
else
  % 3d image
  RoI = SeqLoop.data.RoIZ(:,:,:,1);
  SeqLoop.data.RoIAll = bsxfun(@times, reshape(RoI, [1, size(RoI)]), ...
  [ones(nImages, 1); NaN(size(data.DataAmplitude,1)-nImages, 1)]);
end
data.DataTimeCut = repmat(data.DataTime(1:nImages).', sum(~isnan(RoI(:))), 1);
data.DataAmplitudeCut = reshape(data.DataAmplitude(~isnan(SeqLoop.data.RoIAll)), nImages, []).';

T1starMap = nan(size(RoI));
T1star = nan(sum(~isnan(RoI(:))), 1);
clear fitExpSettings
fitExpSettings.EndOffset = true;
fitExpSettings.RingFilter = false;
fitExpSettings.DoubleExp = false;
fitExpSettings.SingleExpFitType = 1;
fitExpSettings.CorrectPhaseOffset = true;
numDiscarded = 0;
for iPixel = 1:sum(~isnan(RoI(:)))
  time = data.DataTimeCut(iPixel,:);
  amp = data.DataAmplitudeCut(iPixel,:);
  [T, fitExpSettings] = fit_exp(amp, time, fitExpSettings);
  if T.xminSingle(2)/T.xminSingle(1) < -2*SeqLoop.AQSlice(1).RoiFitAmpToOffsetFactorTol ...
      || T.xminSingle(2)/T.xminSingle(1) > -2/SeqLoop.AQSlice(1).RoiFitAmpToOffsetFactorTol
    % Ratio of best-fit amplitude to offset of exponential function deviates
    % from expected value (-2).
    T1star(iPixel) = NaN;
    numDiscarded = numDiscarded + 1;
  else
    T1star(iPixel) = T.tau;
  end
end
fprintf('T1* value for %d out of %d pixels has been discarded because of uncertain inversion in fit.\n', ...
  numDiscarded, numel(T1star));
T1starMap(~isnan(RoI(:))) = T1star;

hf10 = figure(10);
if dim_image < 3
  % 2-d image
  if ~isempty(getappdata(hf10, 'sliceomatic'))
    % previous 3-d results are still open
    clf(hf10);
  end
  ax10(1) = subplot(1,1,1, 'Parent', hf10);
  imagesc(SeqLoop.data.Ticks(1).ReadZ, SeqLoop.data.Ticks(2).PhaseZ, T1starMap.', 'Parent', ax10(1));
  set(ax10(1), 'CLim', [0, SeqLoop.T1Estimated*2]);
  set(ax10(1), 'YDir', 'normal');
  title(ax10(1), {'T1* in s', ...
    sprintf('mean = %.1f ms, STD = %.1f ms', ...
            mean(T1starMap(:), 'omitnan')*1e3, std(T1starMap(:), 'omitnan')*1e3)});
  xlabel(ax10(1), 'Read in m');
  ylabel(ax10(1), 'Phase in m');
  colorbar('peer', ax10(1));
  set(ax10(1), 'DataAspectRatio', [1 1 1]);
  set(hf10, 'Name', 'T1* (Single-Exponential Fit)');
else
  % 3-d image
  hsl = sliceomatic(hf10, permute(squeeze(T1starMap), SeqLoop.data.PermuteOrder), ...
    SeqLoop.data.Ticks(1).ReadZ, SeqLoop.data.Ticks(1).PhaseZ, SeqLoop.data.Ticks(2).PhaseZ);
  title(hsl.hAxes, {'T1* in s', ...
    sprintf('mean = %.1f ms, STD = %.1f ms', ...
            mean(T1starMap(:), 'omitnan')*1e3, std(T1starMap(:), 'omitnan')*1e3)});
  xlabel(hsl.hAxes, [SeqLoop.AQSlice(1).ReadCartesianAxis{1}, 'Read in m']);
  ylabel(hsl.hAxes, [SeqLoop.AQSlice(1).PhaseCartesianAxis{1}, 'Phase(1) in m']);
  zlabel(hsl.hAxes, [SeqLoop.AQSlice(1).PhaseCartesianAxis{2}, 'Phase(2) in m']);
  title(hsl.GetSliderX(), SeqLoop.AQSlice(1).ReadCartesianAxis{1});
  title(hsl.GetSliderY(), SeqLoop.AQSlice(1).PhaseCartesianAxis{1});
  title(hsl.GetSliderZ(), SeqLoop.AQSlice(1).PhaseCartesianAxis{2});
  set(hsl, 'CLim', [0, SeqLoop.T1Estimated*2]);
  set(hf10, 'Name', 'T1* (Single-Exponential Fit)');
  if isa(hsl, 'sliceomatic') && ...
      isempty([hsl.GetAllSlicesPosX(); hsl.GetAllSlicesPosY(); hsl.GetAllSlicesPosZ(); hsl.GetAllIsoValues()])
    % Note: The axes labels don't necessarily correspond to the axes "orientation".
    hsl.AddSliceZ(0);
  end
end


%% Plot data at single pixel
hf9 = figure(9);
set(hf9, 'Name', 'Singular Pixel Data');

iPixel = round(sum(~isnan(RoI(:))/2))+0;
time = data.DataTimeCut(iPixel,:);
amp = data.DataAmplitudeCut(iPixel,:);
fitExpSettingsPixel = fitExpSettings;
fitExpSettingsPixel.hParent = hf9;
fitExpSettingsPixel.DoubleExp = 1;
T = fit_exp(amp, time, fitExpSettingsPixel);


%% -----------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
