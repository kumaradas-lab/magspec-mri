%% Example script to measure a 2d T2 map with a segmented Turbo Spin Echo
%
% Parameters are optimized for an oil sample.
% Do proper shimming with Demo_Auto_ParameterSearch for best results.

LoadSystem;                                               % load system parameters (reset to default: HW Seq AQ TX Grad)

% Define parameters required for the measurement
Seq.T2Estimated = 0.150;                                  % estimated T2 value in seconds
Seq.tEchoTrain  = Seq.T2Estimated*6;                      % duration of echo train in seconds
Seq.LoopsBreak  = max(2,Seq.T2Estimated*5);               % Pause between two loop averages in seconds ([]= fast as possible)

Seq.LoopSaveAllData = 1;                                  % save all data for each image (tEcho)
Seq.LoopPlotAverages = 0;                                 % do not plot average of tEchoes

HW.Grad(1).tRamp = 0.2e-3;                                % increase raise and fall time of gradients
HW.Grad(1).tEC   = 0.1e-3;                                % adjust time delay of gradient pulse

Seq.Loops = 1;                                            % number of loop averages 1...

Seq.tEcho = 7e-3;                                         % echo time in seconds e.g. 5e-3

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 12;                                % number of pixels in read, if nRead>1 nPhase(1)=1
Seq.AQSlice(1).nPhase(1) = 1;                             % number of pixels in phase(1) (3D)
Seq.AQSlice(1).nPhase(2) = 12;                            % number of pixels in phase(2) (2D)
Seq.AQSlice(1).nPhase(3) = 1;                             % number of pixels in phase(3) (CSI)
Seq.AQSlice(1).HzPerPixMin = 1000;                        % bandwidth per pixel in Hz (1/HzPerPixMin= duration of AQ)
Seq.AQSlice(1).sizeRead = 12e-3;                          % size in read direction in meter
Seq.AQSlice(1).sizePhase(1) = 12e-3;                      % size in phase(2) direction in meter
Seq.AQSlice(1).sizePhase(2) = 12e-3;                      % size in phase(2) direction in meter
dim_image = sum([Seq.AQSlice(1).nPhase, Seq.AQSlice(1).nRead] > 1);
if dim_image < 3
  % select slice only for 2-d image
  Seq.AQSlice(1).thickness = 4e-3;                          % slice thickness in meter
  Seq.AQSlice(1).excitationPulse = @Pulse_Sinc_3_Hamming;   % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
else
  % 3-d image
  Seq.AQSlice(1).thickness = Inf;                           % slice thickness in meter
end
Seq.AQSlice(1).inversionPulse = @Pulse_Rect;              % inversion pulse function

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(1) = 1;                            % oversampling phase(1)  1...
Seq.AQSlice(1).PhaseOS(2) = 1;                            % oversampling phase(2)  1...
Seq.AQSlice(1).PhaseOS(3) = 1;                            % oversampling phase(3)  1...
Seq.AQSlice(1).PhaseOS(Seq.AQSlice(1).nPhase==1) = 1;     % reset PhaseOS(x) to 1 if nPhase(x) == 1

% % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).TurboFactor = prod(Seq.AQSlice(1).nPhase(2).*Seq.AQSlice(1).PhaseOS(2))/2;  % number of image k-lines per excitation
Seq.AQSlice(1).TurboBreak = Seq.LoopsBreak;               % break between last echo and next excitation
Seq.tRep = Seq.LoopsBreak;                                % pause between two loop averages in seconds ([]= fast as possible)
Seq.AQSlice(1).oddEvenEchoes = 1;                         % acquire separate images consisting of only odd or even Echoes respectively
Seq.AQSlice(1).phaseCycling = 1;                          % acquire separate images with the inversion pulses rotated by 180 degrees
Seq.AQSlice(1).nImages = ceil(Seq.tEchoTrain/(Seq.tEcho*prod(Seq.AQSlice(1).nPhase(2).*Seq.AQSlice(1).PhaseOS(2)))/(Seq.AQSlice(1).oddEvenEchoes+1));  % number of acquired images

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dim_image < 3
  % 2-d image
  Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xz');  % Top-down image (slice/phase(1), phase(2), read/phase(3))
  % Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zx');  % Bottom-up image (slice/phase(1), phase(2), read/phase(3))
  % Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
  % Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
else
  % 3-d image
  Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zxy');
end

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeqEchoTrain = 1:3;                               % plot sequence wrapped after each echo train, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                         % plot every loop
Seq.AQSlice(1).plotkSpace = 1;                            % plot k-space
Seq.AQSlice(1).plotImage = 1;                             % plot image
Seq.AQSlice(1).plotPhase = 1;                             % plot phase of k-space or image
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;                  % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 2;                        % zero fill resolution factor

Seq.MaxGradAmpSlice = 0.02;                               % limit slice gradient strength
Seq.CorrectSliceRephase = 0*double(Seq.AQSlice(1).thickness<=0.015); % correct SliceGradTimeIntegralRephaseOffset
Seq.CorrectPhaseRephase = 0;                              % correct PhaseGradTimeIntegralRephaseOffset
Seq.CorrectReadRephase = 0;                               % correct ReadGradTimeIntegralOffset
Seq.AQSlice(1).SpoilFactor = [0, 0.00, Inf];              % spoiler resolution factor (slice or phase 1, phase 2, read or phase 3)
Seq.AQSlice(1).SpoilLengthFactor = 1;                     % spoil length factor

Seq.SteadyState_PreShots = 0;
Seq.SteadyState_PreShots180=0;

HW.FindFrequencyPause = Seq.T2Estimated*5;                % increase pause after frequency search
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 1);       % find magnet frequency
Seq.Find_Frequency_interval = 1000;                       % run Find_Frequency every X seconds or after every loop (if measurement duration > X sec)

[SeqLoop, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);  % run measurement


if 0
  %% Set new Zero Fill parameters
  SeqLoop.AQSlice(1).ZeroFillWindowSize = 1.11;               % zero fill window size (high k-space values are damped by a cos^2 law)
  SeqLoop.AQSlice(1).ZeroFillFactor = 2;                      % zero fill resolution factor
  SeqLoop.data.RoI = [];
  SeqLoop.data = get_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
end


%% Collect and process data from SeqLoop
clear data
% get images
data.ImageZ = SeqLoop.data.ImageZ;
% phase correction with first image
data.ImageZ = bsxfun(@times, data.ImageZ, exp(-1i*(0+angle(data.ImageZ(:,:,:,:,1)))));
% subtract last image to reduce artifacts from stimulated Echoes
data.ImageZ = bsxfun(@minus, data.ImageZ, data.ImageZ(:,:,:,:,end));
% data.ImageZ(:,:,:,:,1:end-1) = bsxfun(@minus, data.ImageZ(:,:,:,:,1:end-1), data.ImageZ(:,:,:,:,end));
% get k-space of all images
data.kSpaceOsRaw = SeqLoop.data.kSpaceOsRaw;
% phase correction with first k-space
data.kSpaceOsRaw = bsxfun(@times, data.kSpaceOsRaw, exp(-1i*(0+angle(data.kSpaceOsRaw(:,:,:,:,1)))));
% subtract last k-space to reduce artifacts from stimulated Echoes
data.kSpaceOsRaw = bsxfun(@minus, data.kSpaceOsRaw, data.kSpaceOsRaw(:,:,:,:,end));
% data.kSpaceOsRaw(:,:,:,:,1:end-1)=bsxfun(@minus, data.kSpaceOsRaw(:,:,:,:,1:end-1),data.kSpaceOsRaw(:,:,:,:,end));

amax = max(abs(data.ImageZ(:)));

% permute dimension with tEchoes to front in amplitude
data.DataAmplitude = permute(real(data.ImageZ(:,:,:,:,1:end-1)), [5,1,2,3,4])/amax;
data.DataTime = mean(SeqLoop.data.tImageZ(1:end-1,:), 2);

% sliceomatic(1, squeeze(abs(data.ImageZ))) % Plot T2 image stack
% Plot some images of the T2 image stack
minmax = [-max(abs(data.ImageZ(:))), max(abs(data.ImageZ(:)))] / amax;
amaxk = max(abs(data.kSpaceOsRaw(:)));
minmaxk = [-max(abs(data.kSpaceOsRaw(:))), max(abs(data.kSpaceOsRaw(:)))] / amaxk;
numAll = ceil([1, SeqLoop.AQSlice(1).nImages/81, SeqLoop.AQSlice(1).nImages/27, SeqLoop.AQSlice(1).nImages/9, SeqLoop.AQSlice(1).nImages/3, SeqLoop.AQSlice(1).nImages-1]);
if sum(numAll==1) == 2
  numAll = numAll(2:end);
elseif sum(numAll==1) > 2
  if SeqLoop.AQSlice(1).nImages < 8
    numAll = 1:SeqLoop.AQSlice(1).nImages-1;
  elseif SeqLoop.AQSlice(1).nImages < 16
    numAll = [1,2,3,6:3:SeqLoop.AQSlice(1).nImages-1];
  elseif SeqLoop.AQSlice(1).nImages < 24
    numAll = [1,2,4,7:4:SeqLoop.AQSlice(1).nImages-1];
  elseif SeqLoop.AQSlice(1).nImages < 8*5
    numAll = [1,2,4,9:6:SeqLoop.AQSlice(1).nImages-1];
  elseif SeqLoop.AQSlice(1).nImages < 8*6
    numAll=[1,2,5,12:7:SeqLoop.AQSlice(1).nImages-1];
  else % if SeqLoop.AQSlice(1).nImages < 57
    numAll = [1,2,6,13,20:9:SeqLoop.AQSlice(1).nImages-1];
  end
  numAll(end) = SeqLoop.AQSlice(1).nImages-1;
end


iAxes = 0;
dim_image = sum([SeqLoop.AQSlice(1).nPhase, SeqLoop.AQSlice(1).nRead] > 1);
if dim_image < 3
  % 2-d image
  if Seq.AQSlice(1).plotkSpace
    k = 2;
  else
    k = 1;
  end
  hf = figure(2);
  clf(hf, 'reset');

  for num = numAll
    iAxes = iAxes+1;
    minmax2 = [-max(max(max(max(max(abs(data.ImageZ(:,:,:,:,num)))))))-eps(amax), ...
      max(max(max(max(max(abs(data.ImageZ(:,:,:,:,num)))))))+eps(amax)]/amax;
    % image amplitude
    hax = subplot(numel(numAll),5*k,1+(iAxes-1)*5*k);
    imagesc(squeeze(abs(data.ImageZ(:,:,:,:,num))).'/amax, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmax2);
    title(hax, ['abs ' num2str(1/minmax2(2),'%.1f')]);
    ylabel(hax, ['TE' num2str(num) '=' num2str(data.DataTime(num)),' s']);
    % image real
    hax = subplot(numel(numAll),5*k,2+(iAxes-1)*5*k, 'Parent', hf);
    imagesc(squeeze(real(data.ImageZ(:,:,:,:,num))).'/amax, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmax2);
    title(hax, ['real ' num2str(1/minmax2(2),'%.1f')]);
    % image imag
    hax = subplot(numel(numAll),5*k,3+(iAxes-1)*5*k, 'Parent', hf);
    imagesc(squeeze(imag(data.ImageZ(:,:,:,:,num))).'/amax, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmax2);
    title(hax, ['imag ' num2str(1/minmax2(2),'%.1f')]);
    % image phase
    hax = subplot(numel(numAll),5*k,4+(iAxes-1)*5*k, 'Parent', hf);
    imagesc(squeeze(angle(data.ImageZ(:,:,:,:,num))).', 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', [-pi,pi]);
    title(hax, 'phase');
    % image amplitude with color limits of first row
    hax = subplot(numel(numAll),5*k,5+(iAxes-1)*5*k, 'Parent', hf);
    imagesc(squeeze(abs(data.ImageZ(:,:,:,:,num))).'/amax, 'Parent', hax);
    set(hax,'YDir','normal', 'CLim', minmax);
    title(hax, ['abs ' num2str(1/minmax2(2),'%.1f') ', const. CLim']);
    if k == 2
      minmaxk2 = [-max(max(max(max(max(abs(data.kSpaceOsRaw(:,:,:,:,num)))))))-eps(amaxk), ...
        max(max(max(max(max(abs(data.kSpaceOsRaw(:,:,:,:,num)))))))+eps(amaxk)]/amaxk;
      % kSpace amplitude
      hax = subplot(numel(numAll),5*k,1+(iAxes-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(abs(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', minmaxk2);
      title(hax, ['abs ' num2str(1/minmaxk2(2),'%.1f')]);
      % kSpace real
      hax = subplot(numel(numAll),5*k,2+(iAxes-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(real(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', minmaxk2);
      title(hax, ['real ' num2str(1/minmaxk2(2),'%.1f')]);

      hax = subplot(numel(numAll),5*k,3+(iAxes-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(imag(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', minmaxk2);
      title(hax, ['imag ' num2str(1/minmaxk2(2),'%.1f')]);

      hax = subplot(numel(numAll),5*k,4+(iAxes-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(angle(data.kSpaceOsRaw(:,:,:,:,num))).', 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', [-pi,pi]);
      title(hax, 'phase');

      hax = subplot(numel(numAll),5*k,5+(iAxes-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(abs(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
      set(hax,'YDir','normal', 'CLim', minmaxk);
      title(hax, 'abs 1.0, const. CLim');
    end

  end

  drawnow();
end

if isemptyfield(SeqLoop.AQSlice(1), 'RoiRelativeValue')
  SeqLoop.AQSlice(1).RoiRelativeValue = 1/5;
end
if isemptyfield(SeqLoop.AQSlice(1), 'RoiCutOffPercentile')
  SeqLoop.AQSlice(1).RoiCutOffPercentile = .95;
end
imageSorted = sort(abs(SeqLoop.data.ImageZ(:)));
SeqLoop.data.RoICutOff = imageSorted(round(numel(SeqLoop.data.ImageZ)*SeqLoop.AQSlice(1).RoiCutOffPercentile));
clear imageSorted
SeqLoop.data.RoI = ones(size(squeeze(SeqLoop.data.ImageZ)));
SeqLoop.data.RoI(abs(squeeze(SeqLoop.data.ImageZ))<=SeqLoop.data.RoICutOff.*SeqLoop.AQSlice(1).RoiRelativeValue) = 0;
SeqLoop.data.RoI = convn(SeqLoop.data.RoI, ones(3,3), 'same');
SeqLoop.data.RoI(SeqLoop.data.RoI<2*3) = NaN;
SeqLoop.data.RoI(~isnan(SeqLoop.data.RoI)) = 1;


%% Single-exponential fit of amplitude in each pixel
nT2Estimated = 3;

n = min(size(data.DataAmplitude,1), ...
  ceil(nT2Estimated*SeqLoop.T2Estimated/diff(data.DataTime(1:2))));
if dim_image < 3
  % 2-d image
  RoI = SeqLoop.data.RoI(:,:,1);
  SeqLoop.data.RoIAll = bsxfun(@times, ...
    reshape(RoI, 1, size(RoI, 1), 1, size(RoI, 2)), ...
    [ones(n, 1); NaN(size(data.DataAmplitude,1)-n, 1)]);
else
  % 3-d image
  RoI = SeqLoop.data.RoI(:,:,:,1);
  SeqLoop.data.RoIAll = bsxfun(@times, ...
    reshape(RoI, [1, size(RoI)]), ...
    [ones(n, 1); NaN(size(data.DataAmplitude,1)-n, 1)]);
end

data.DataTimeCut = repmat(data.DataTime(1:n).', sum(~isnan(RoI(:))), 1);
data.DataAmplitudeCut = abs(reshape(data.DataAmplitude(~isnan(SeqLoop.data.RoIAll)), n, [])).';


T2map = NaN(size(RoI));
T2 = NaN(sum(~isnan(RoI(:))), 1);
clear fitExpSettings
fitExpSettings.EndOffset = false;
fitExpSettings.DoubleExp = false;
fitExpSettings.SingleExpFitType = 2;
for iPixel = 1:sum(~isnan(RoI(:)))
  time = data.DataTimeCut(iPixel,:);
  amp = data.DataAmplitudeCut(iPixel,:);
  [T, fitExpSettings] = fit_exp(amp, time, fitExpSettings);
  T2(iPixel) = T.tau;
end
T2map(~isnan(RoI(:))) = T2;


hf10 = figure(10);
if dim_image < 3
  % 2-d image
  if ~isempty(getappdata(hf10, 'sliceomatic'))
    % previous 3d results are still open
    clf(hf10);
  end
  ax10(1) = subplot(1,1,1, 'Parent', hf10);
  imagesc(SeqLoop.data.Ticks(1).ReadZ, SeqLoop.data.Ticks(2).PhaseZ, T2map.', 'Parent', ax10(1));
  set(ax10(1), 'CLim', [0,SeqLoop.T2Estimated*2]);
  set(ax10(1), 'YDir', 'normal');
  title(ax10(1), ['T2 in s, mean = ' num2str(mean(T2)) ' s, STD = ' num2str(std(T2)) ' s']);
  xlabel(ax10(1), 'Read in m');
  ylabel(ax10(1), 'Phase in m');
  colorbar('peer', ax10(1));
  set(ax10(1), 'DataAspectRatio', [1 1 1]);
  set(hf10, 'Name', 'T2 (Single-Exponential Fit)');
else
  % 3-d image
  hsl = sliceomatic(hf10, permute(squeeze(T2map), SeqLoop.data.PermuteOrder), SeqLoop.data.Ticks(1).ReadZ, SeqLoop.data.Ticks(1).PhaseZ, SeqLoop.data.Ticks(2).PhaseZ);
  title(hsl.hAxes, {'T2 in s', ['mean = ' num2str(mean(T2)) ' s, STD = ' num2str(std(T2)) ' s']});
  xlabel(hsl.hAxes, [SeqLoop.AQSlice(1).ReadCartesianAxis{1}, 'Read in m']);
  ylabel(hsl.hAxes, [SeqLoop.AQSlice(1).PhaseCartesianAxis{1}, 'Phase(1) in m']);
  zlabel(hsl.hAxes, [SeqLoop.AQSlice(1).PhaseCartesianAxis{2}, 'Phase(2) in m']);
  title(hsl.GetSliderX(), SeqLoop.AQSlice(1).ReadCartesianAxis{1});
  title(hsl.GetSliderY(), SeqLoop.AQSlice(1).PhaseCartesianAxis{1});
  title(hsl.GetSliderZ(), SeqLoop.AQSlice(1).PhaseCartesianAxis{2});
  set(hf10, 'Name', 'T2 (Single-Exponential Fit)');
  if isa(hsl, 'sliceomatic') && ...
      isempty([hsl.GetAllSlicesPosX(); hsl.GetAllSlicesPosY(); hsl.GetAllSlicesPosZ(); hsl.GetAllIsoValues()])
    % Note: The axes labels don't necessarily correspond to the axes "orientation".
    hsl.AddSliceX(0);
  end
end


%% Plot data at single pixel
hf9 = figure(9);
set(hf9, 'Name', 'Single pixel data');

iPixel = round(sum(~isnan(RoI(:))/2))+0;
time = data.DataTimeCut(iPixel,:);
amp = data.DataAmplitudeCut(iPixel,:);
fitExpSettingsPixel = fitExpSettings;
fitExpSettingsPixel.hParent = hf9;
fitExpSettingsPixel.DoubleExp = 1;
T = fit_exp(amp, time, fitExpSettingsPixel);


%% -----------------------------------------------------------------------------
% (C) Copyright 2019-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
