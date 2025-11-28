%% Example script to measure a 2d or 3d T2 map with Spin Echo

LoadSystem;                                                 % load system parameters (reset to default: HW Seq AQ TX Grad)
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 1);         % find magnet frequency

Seq.T2Estimated = 0.12;                                     % estimated T2 value in seconds, e.g. 120ms for Oil
Seq.T1Estimated = Seq.T2Estimated;                          % estimated T1 value in seconds, e.g. 120ms for Oil
Seq.tEchoTrain  = Seq.T2Estimated*6;                        % duration of Echo train
Seq.LoopsBreak  = Seq.T1Estimated*3;                        % pause between two loop averages in seconds ([]= fast as possible)

Seq.LoopSaveAllData = 1;
Seq.LoopPlotAverages = 0;

Seq.Loops = 1;                                              % number of loop averages    1...
% Seq.Find_Frequency_interval = 100;                          % time between two frequency searches
Seq.plotSeq = 1:3;                                          % plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
% Seq.plotSeqStart = 10;                                      % plot sequence start with tRep x
% Seq.plotSeqEnd = 11;                                        % plot sequence stop with tRep x
Seq.plot3D = 0;                                             % plot iLaplace T2 fit

Seq.tEcho = 2.5e-3;                                         % echo time in seconds e.g. 5e-3
Seq.tRep = Seq.LoopsBreak;                                  % repetition time in seconds; only used if Seq.TurboBreak=[], auto Seq.tRep=[] and Seq.AQSlice(1).TurboBreak= x sec;

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 12;                                  % number of pixels in read, if nRead>1 nPhase(1)=1
Seq.AQSlice(1).nPhase(1) = 1;                               % number of pixels in phase(1)
Seq.AQSlice(1).nPhase(2) = 12;                              % number of pixels in phase(2)
Seq.AQSlice(1).HzPerPixMin = 2000;                          % bandwith per pixel in Hz (1/HzPerPixMin = duration of AQ)
Seq.AQSlice(1).sizeRead = 0.012;                            % image size in read in meter (for CSI set to 1e12)
Seq.AQSlice(1).sizePhase(1) = 0.012;                        % image size in phase(1) in meter
Seq.AQSlice(1).sizePhase(2) = 0.012;                        % image size in phase(2) in meter
dim_image = sum([Seq.AQSlice(1).nPhase, Seq.AQSlice(1).nRead] > 1);
if dim_image < 3
  % select slice only for 2d image
  Seq.AQSlice(1).thickness = 0.003;                         % image thickness in slice direction  used for 2D and 3D! ([] for no Slice) in meter
  Seq.AQSlice(1).thicknessInversion = Seq.AQSlice(1).thickness;  % image thickness in slice direction  used for 2D and 3D! ([] for no Slice) in meter
end
Seq.AQSlice(1).excitationPulse = @Pulse_Rect;               % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.AQSlice(1).inversionPulse = @Pulse_Rect;                % inversion pulse function

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.AQSlice(1).ReadOS = 16;                                 % oversampling read ([] for automatic, recommended >=2) 1...
Seq.AQSlice(1).PhaseOS(1) = 2;                              % oversampling phase(1)  1...
Seq.AQSlice(1).PhaseOS(2) = 2;                              % oversampling phase(2)  1...
Seq.AQSlice(1).PhaseOS(3) = 1;                              % oversampling phase(3)  1... (set to 1 if nPhase(3)=1;
Seq.AQSlice(1).PhaseOS(Seq.AQSlice(1).nPhase==1) = 1;       % reset PhaseOS(x) to 1 if nPhase(x) == 1

% % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).oddEvenEchoes = 0;                           % acquire separate images consisting of only odd or even Echoes respectively
Seq.AQSlice(1).phaseCycling = 1;                            % acquire separate images with the inversion pulses rotated by 180 degrees (supress Echoes of inversion pulses)
Seq.AQSlice(1).nImages = round(Seq.tEchoTrain/Seq.tEcho/(Seq.AQSlice(1).oddEvenEchoes+1));   % number of acquired images
Seq.AQSlice(1).TurboFactor = 1;                             % number of image k-lines per excitation
Seq.AQSlice(1).TurboBreak = Seq.tRep;                       % break between last echo and next excitation

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.LoopPlot = 1;                                           % plot every loop
Seq.AQSlice(1).plotkSpace = 1;                              % plot k-space
Seq.AQSlice(1).plotImage = 1;                               % plot image
Seq.AQSlice(1).plotPhase = 1;                               % plot phase of k-space and/or image
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;                    % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 2;                          % zero fill resolution factor

% % Orientation in Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dim_image < 3
  % 2d image
  Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zx');
else
  % 3d image
  Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zxy');
end
% Seq.AQSlice(1).alfa = 0.0*pi;                               % first rotation around x axis in RAD
% Seq.AQSlice(1).phi = 0.0*pi;                                % second rotation around y axis in RAD
% Seq.AQSlice(1).theta = -0.5*pi;                             % third rotation around z axis in RAD
Seq.AQSlice(1).Center2OriginImage = [0.00,0.000,0.000];     % vector from center of the image to the origin in Image coordinate system

% % Some corrections  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.AQSlice(1).SpoilLengthFactor = 1;                       % spoil length factor
% Seq.AQSlice(1).sizePhaseSpoil = [1,1,1]*0.0005;%1e12;       % spoil size
Seq.AQSlice(1).SpoilFactor = [4,0,0];                       % spoiler resolution factor (slice or phase 1, phase 2, read or phase 3)
Seq.MaxGradAmpSlice = 0.05;                                 % limit slice gradient strength
Seq.SteadyState_PreShots90 = 0;                             % number of echo trains before first AQ
Seq.SteadyState_PreShots180 = 0;                            % number of echoes before first AQ

Seq.PhaseCorrectionWithFirstImage = 1;                      % phase correction with first image and k-space
Seq.SubtractLastImage = 0;                                  % subtract last image to reduce artifacts from stimulated Echoes


[SeqLoop, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);

if 0
  %% Set new Zero Fill parameter
  SeqLoop.AQSlice(1).ZeroFillWindowSize = 1.11;              % zero fill window size (high k-space values are damped by a cos^2 law)
  SeqLoop.AQSlice(1).ZeroFillFactor = 2;                     % zero fill resolution factor
  SeqLoop.data = get_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
end


%% Collect and process data from SeqLoop
clear data
data.ImageZ = SeqLoop.data.ImageZ; % get images
data.kSpaceOsRaw = SeqLoop.data.kSpaceOsRaw; % get k-space of all images

if SeqLoop.PhaseCorrectionWithFirstImage
  data.ImageZFirstPhaseCor = exp(-1i*(0+angle(data.ImageZ(:,:,:,:,1))));
  data.kSpaceOsRawFirstPhaseCor = exp(-1i*(0+angle(data.kSpaceOsRaw(:,:,:,:,1))));
  % phase correction with first image and k-space
  data.ImageZ = bsxfun(@times, data.ImageZ, exp(-1i*(0+angle(data.ImageZ(:,:,:,:,1)))));
  data.kSpaceOsRaw = bsxfun(@times, data.kSpaceOsRaw, exp(-1i*(0+angle(data.kSpaceOsRaw(:,:,:,:,1)))));
else
  data.ImageZFirstPhaseCor = ones(size(data.ImageZFirstPhaseCor));
  data.kSpaceOsRawFirstPhaseCor = ones(size(data.kSpaceOsRawFirstPhaseCor));
end

amax = max(abs(data.ImageZ(:)));
if SeqLoop.SubtractLastImage
  % subtract last image and k-space to reduce artifacts from stimulated Echoes
  data.ImageZ = bsxfun(@minus, data.ImageZ, data.ImageZ(:,:,:,:,end));  % subtract last image
  data.kSpaceOsRaw = bsxfun(@minus, data.kSpaceOsRaw, data.kSpaceOsRaw(:,:,:,:,end)); % subtract last kSpace
  % permute dimension with tEchoes to front in amplitude
  data.DataAmplitude = permute(data.ImageZ(:,:,:,:,1:end-1), [5,1,2,3,4])/amax;
  data.DataTime = mean(SeqLoop.data.tImageZ(1:end-1,:), 2);
  SeqLoop.nLastImage = SeqLoop.AQSlice(1).nImages-1;
else
  % permute dimension with tEchoes to front in amplitude
  data.DataAmplitude = permute(data.ImageZ(:,:,:,:,1:end), [5,1,2,3,4])/amax;
  data.DataTime = mean(SeqLoop.data.tImageZ(1:end,:), 2);
  SeqLoop.nLastImage = SeqLoop.AQSlice(1).nImages;
end

% sliceomatic(1, squeeze(abs(data.ImageZAll))) % Plot T2 image stack
% Plot some images of the T2 image stack
iAxes = 0;
minmax = [-max(abs(data.ImageZ(:))),max(abs(data.ImageZ(:)))]/amax;
amaxk = max(abs(data.kSpaceOsRaw(:)));
minmaxk = [-max(abs(data.kSpaceOsRaw(:))), max(abs(data.kSpaceOsRaw(:)))]/amaxk;
numAll = ceil([1, SeqLoop.AQSlice(1).nImages/81, SeqLoop.AQSlice(1).nImages/27, SeqLoop.AQSlice(1).nImages/9, SeqLoop.AQSlice(1).nImages/3, SeqLoop.nLastImage]);
if sum(numAll==1) == 2
  numAll = numAll(2:end);
elseif sum(numAll==1) > 2
  if SeqLoop.AQSlice(1).nImages < 8
    numAll = 1:SeqLoop.nLastImage;
  elseif SeqLoop.AQSlice(1).nImages < 16
    numAll = [1,2,3,6:3:SeqLoop.nLastImage];
  elseif SeqLoop.AQSlice(1).nImages < 24
    numAll = [1,2,4,7:4:SeqLoop.nLastImage];
  elseif SeqLoop.AQSlice(1).nImages < 8*5
    numAll = [1,2,4,9:6:SeqLoop.nLastImage];
  elseif SeqLoop.AQSlice(1).nImages < 8*6
    numAll=[1,2,5,12:7:SeqLoop.nLastImage];
  else % if SeqLoop.AQSlice(1).nImages < 57
    numAll = [1,2,6,13,20:9:SeqLoop.nLastImage];
  end
  numAll(end) = SeqLoop.AQSlice(1).nImages-1;
end
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
    title(hax, ['abs ' num2str(1/minmax2(2),'%.1f')]);
    ylabel(hax, ['TE' num2str(num) '=' num2str(SeqLoop.data.tImageZ(num)),' s']);
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
    if k==2
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
      title(hax, 'abs 1.00');
    end

  end
  drawnow();
end

SeqLoop.AQSlice = set_EmptyField(SeqLoop.AQSlice(1), 'RoiRelativeValue', 1/5);
SeqLoop.AQSlice = set_EmptyField(SeqLoop.AQSlice(1), 'RoiCutOffPercentile', .95);
imageSorted = sort(abs(SeqLoop.data.ImageZ(:)));
SeqLoop.data.RoICutOff = imageSorted(round(numel(SeqLoop.data.ImageZ)*SeqLoop.AQSlice(1).RoiCutOffPercentile));
clear imageSorted
SeqLoop.data.RoI = ones(size(squeeze(SeqLoop.data.ImageZ)));
SeqLoop.data.RoI(abs(squeeze(SeqLoop.data.ImageZ))<=SeqLoop.data.RoICutOff.*SeqLoop.AQSlice(1).RoiRelativeValue) = 0;
SeqLoop.data.RoI = convn(SeqLoop.data.RoI, ones(3,3), 'same');
SeqLoop.data.RoI(SeqLoop.data.RoI<2*3) = nan;
SeqLoop.data.RoI(~isnan(SeqLoop.data.RoI)) = 1;

%% Single-exponential fit of amplitude in each pixel
nT2Estimated = 6;

n = min(size(data.DataAmplitude,1), ceil(nT2Estimated*SeqLoop.T2Estimated/diff(data.DataTime(1:2))));
if dim_image < 3
  % 2d image
  RoI = SeqLoop.data.RoI(:,:,1);
  SeqLoop.data.RoIAll = bsxfun(@times, reshape(RoI, 1, size(RoI,1), 1, size(RoI,2)), [ones(n,1);nan(size(data.DataAmplitude,1)-n,1)]);
else
  % 3d image
  RoI = SeqLoop.data.RoI(:,:,:,1);
  SeqLoop.data.RoIAll = bsxfun(@times, reshape(RoI, [1, size(RoI)]), [ones(n,1);nan(size(data.DataAmplitude,1)-n,1)]);
end
data.DataTimeCut = repmat(data.DataTime(1:n).', sum(~isnan(RoI(:))), 1);
data.DataAmplitudeCut = reshape(data.DataAmplitude(~isnan(SeqLoop.data.RoIAll)), n, []).';

T2map = nan(size(RoI));
T2 = nan(sum(~isnan(RoI(:))), 1);
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
  % 2d image
  if ~isempty(getappdata(hf10, 'sliceomatic'))
    % previous 3d results are still open
    clf(hf10);
  end
  ax10(1) = subplot(1,1,1, 'Parent', hf10);
  imagesc(SeqLoop.data.Ticks(1).ReadZ, SeqLoop.data.Ticks(2).PhaseZ, T2map.', 'Parent', ax10(1));
  set(ax10(1), 'CLim', [0,SeqLoop.T2Estimated*2]);
  set(ax10(1), 'YDir', 'normal');
  title(ax10(1), {'T2 in s', ['mean = ' num2str(mean(T2)) ' s, STD = ' num2str(std(T2)) ' s']});
  xlabel(ax10(1), 'Read in m');
  ylabel(ax10(1), 'Phase in m');
  colorbar('peer', ax10(1));
  set(ax10(1), 'DataAspectRatio', [1 1 1]);
  set(hf10, 'Name', 'T2 (Single-Exponential Fit)');
else
  % 3d image
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

%% Plot data at singular pixel
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
% (C) Copyright 2016-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
