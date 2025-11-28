%% Gradient Echo Train 2D (Flash 2D)
% This example acquires multiple 2D images from gradient echo trains with Ernst
% angle excitation.
% Subsequently, an exponential decay curve is fit to each pixel/voxel in the
% image stack. The result is a combination of the T2star of the sample and the
% B0 inhomogeneity.

%%
LoadSystem;                                         % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.Loops = 1;                                      % number of loop averages 1...

Seq.T1 = 100e-3;                                    % T1 of sample for excitation with Ernst angle: acosd(exp(-Seq.tRep/Seq.T1))
Seq.tEcho = 4e-3;                                   % echo time in seconds e.g. 4e-3
Seq.tRep = 200e-3;                                  % repetition time in seconds (default is Seq.tEcho*2)

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 32;                          % number of pixels in read direction
Seq.AQSlice(1).nPhase(2) = 32;                      % number of pixels in phase direction
Seq.AQSlice(1).HzPerPixMin = 0;                     % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ window, 0: longest possible)
Seq.AQSlice(1).sizeRead = 0.020;                    % size in read direction in meter
Seq.AQSlice(1).sizePhase(2) = 0.020;                % size in phase(2) direction in meter
Seq.AQSlice(1).thickness = 0.003;                   % slice thickness in meter
Seq.AQSlice(1).excitationPulse = @Pulse_RaisedCos;  % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.MaxGradAmpSlice = 0.04;                         % maximum gradient amplitude for slice selection in T/m (reduce to lessen impact of eddy currents)

% % EPI settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nImages = 12;                        % number of images (k-lines) acquired in one gradient echo train
Seq.AQSlice(1).EPIGradient = 1;                     % move encoding gradient to allow them to cancel each other
Seq.AQSlice(1).EPIReadBipolar = 1;                  % Use bipolar read gradients in gradient echo train

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(2) = 2;                      % oversampling phase(2)  1...

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xz');  % Top-down image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zx');  % Bottom-up image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeqAQ = 1:3;                                % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                   % plot every loop
Seq.AQSlice(1).plotkSpace = 1;                      % plot k-space
Seq.AQSlice(1).plotImage = 1;                       % plot image
Seq.AQSlice(1).plotPhase = 1;                       % plot phase of k-space and/or image
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;            % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 4;                  % zero fill resolution factor

Seq.CorrectSliceRephase = 0;                        % correct SliceGradTimeIntegralOffset
Seq.CorrectReadRephase = 0;                         % correct ReadGradTimeIntegralOffset
Seq.CorrectPhase = 1;                               % acquire FID directly after FID to track frequency
Seq.CorrectPhaseDuration = 1.5e-3;                  % duration of frequency tracking windows in s
Seq.CorrectPhase_maxPhaseOffset = 2*pi/10;          % maximum tolerated phase offset in rad due to standard error of frequency tracking

% % Post processing settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.PhaseCorrectionWithFirstImage = 1;              % phase correction with first image and k-space
Seq.SubtractLastImage = 0;                          % subtract last image to reduce artifacts from stimulated Echoes


[SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);

if 0
  %% Test zero fill and k-space filter window
  % SeqLoop.AQSlice(1).ZeroFillWindowSize = 1.4;      % zero fill window size (high k-space values are damped by a cos^2 law)
  % SeqLoop.AQSlice(1).ZeroFillFactor = 8;            % zero fill resolution factor
  % SeqLoop.data.RoI = [];                            % RoI has to be reset in order to recalculate
  % SeqLoop.AQSlice(1).RoiCutOffPercentile = [];      % percentile that is the reference for the ROI selection
  % SeqLoop.AQSlice(1).RoiRelativeValue = 0;          % relative value of that percentile
  SeqLoop.AQSlice(1).iSlice = 3;                    % number of the image in the echo train to show

  [SeqLoop.data] = get_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
  [SeqLoop.data, SeqLoop.AQSlice] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
end

%% Collect and process data from SeqLoop
clear data
data.ImageZ = SeqLoop.data.ImageZ;  % get images
data.kSpaceOsRaw = SeqLoop.data.kSpaceOsRaw;  % get k-space of all images

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
  data.kSpaceOsRaw = bsxfun(@minus, data.kSpaceOsRaw, data.kSpaceOsRaw(:,:,:,:,end));  % subtract last kSpace
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

% sliceomatic(1, squeeze(abs(data.ImageZAll)));  % Plot T2star image stack
% Plot some images of the T2star image stack
iAxes = 0;
hf = figure(2);
clf(hf, 'reset')
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
    numAll = [1, 2, 3, 6:3:SeqLoop.nLastImage];
  elseif SeqLoop.AQSlice(1).nImages < 24
    numAll = [1, 2, 4, 7:4:SeqLoop.nLastImage];
  elseif SeqLoop.AQSlice(1).nImages < 8*5
    numAll = [1, 2, 4, 9:6:SeqLoop.nLastImage];
  elseif SeqLoop.AQSlice(1).nImages < 8*6
    numAll = [1, 2, 5, 12:7:SeqLoop.nLastImage];
  else % if SeqLoop.AQSlice(1).nImages < 57
    numAll = [1, 2, 6, 13, 20:9:SeqLoop.nLastImage];
  end
  numAll(end) = SeqLoop.AQSlice(1).nImages-1;
end
if Seq.AQSlice(1).plotkSpace, k=2; else k=1; end
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
  set(hax, 'YDir', 'normal', 'CLim', [-pi, pi]);
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
    set(hax, 'YDir', 'normal', 'CLim', [-pi, pi]);
    title(hax, 'phase');

    hax = subplot(numel(numAll),5*k,5+(iAxes-1)*5*k+5, 'Parent', hf);
    imagesc(squeeze(abs(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
    set(hax,'YDir','normal', 'CLim', minmaxk);
    title(hax, 'abs 1.00');
  end

end
drawnow();

SeqLoop.AQSlice = set_EmptyField(SeqLoop.AQSlice(1), 'RoiRelativeValue', 1/5);
SeqLoop.AQSlice = set_EmptyField(SeqLoop.AQSlice(1), 'RoiCutOffPercentile', .95);
imageSorted = sort(abs(SeqLoop.data.ImageZ(:)));
SeqLoop.data.RoICutOff = imageSorted(round(numel(SeqLoop.data.ImageZ)*SeqLoop.AQSlice(1).RoiCutOffPercentile));
clear imageSorted
SeqLoop.data.RoI = ones(size(squeeze(SeqLoop.data.ImageZ)));
SeqLoop.data.RoI(abs(squeeze(SeqLoop.data.ImageZ))<=SeqLoop.data.RoICutOff.*SeqLoop.AQSlice(1).RoiRelativeValue) = 0;
SeqLoop.data.RoI = convn(SeqLoop.data.RoI, ones(3,3), 'same');
SeqLoop.data.RoI(SeqLoop.data.RoI<2*3) = NaN;
SeqLoop.data.RoI(~isnan(SeqLoop.data.RoI)) = 1;

%% Single-exponential fit of amplitude in each pixel
nT2Estimated = SeqLoop.AQSlice(1).nImages;

n = min(size(data.DataAmplitude,1), ceil(SeqLoop.tEcho(end)/diff(data.DataTime(1:2))));
RoI = SeqLoop.data.RoI(:,:,1);
SeqLoop.data.RoIAll = bsxfun(@times, reshape(RoI, 1, size(RoI,1), 1, size(RoI,2)), ...
  [ones(n,1);nan(size(data.DataAmplitude,1)-n,1)]);

data.DataTimeCut = repmat(data.DataTime(1:n).', sum(~isnan(RoI(:))), 1);
data.DataAmplitudeCut = reshape(data.DataAmplitude(~isnan(SeqLoop.data.RoIAll)), n, []).';

T2map = nan(size(RoI));
T2star = nan(sum(~isnan(RoI(:))), 1);
clear fitExpSettings
fitExpSettings.EndOffset = false;
fitExpSettings.RingFilter = false;
fitExpSettings.CorrectFrequencyOffset = true;
fitExpSettings.CorrectPhaseOffset = true;
fitExpSettings.DoubleExp = false;
fitExpSettings.SingleExpFitType = 2;
for iPixel = 1:sum(~isnan(RoI(:)))
  time = data.DataTimeCut(iPixel,:);
  amp = data.DataAmplitudeCut(iPixel,:);
  [T, fitExpSettings] = fit_exp(amp*amax, time, fitExpSettings);
  T2star(iPixel) = T.tau;
end
T2map(~isnan(RoI(:))) = T2star;

hf10 = figure(10);
ax10(1) = subplot(1,1,1, 'Parent', hf10);
imagesc(SeqLoop.data.Ticks(1).ReadZ, SeqLoop.data.Ticks(2).PhaseZ, T2map.', 'Parent', ax10(1));
% minimum value for color axis limits
maxCLim = max([mean(T2star(:), 'omitnan')*2, SeqLoop.tEcho(end)/4]);
set(ax10(1), 'CLim', [0, maxCLim]);
set(ax10(1), 'YDir', 'normal');
title(ax10(1), {'T2* in s', ['mean = ' num2str(mean(T2star)) ' s, \sigma = ' num2str(std(T2star)) ' s']});
xlabel(ax10(1), 'Read in m');
ylabel(ax10(1), 'Phase in m');
colorbar('peer', ax10(1));
set(ax10(1), 'DataAspectRatio', [1, 1, 1]);
set(hf10, 'Name', 'T2* (Single-Exponential Fit)');


%% Plot data at singular pixel
hf9 = figure(9);
set(hf9, 'Name', 'Singular Pixel Data');

iPixel = round(sum(~isnan(RoI(:))/2))+17;
time = data.DataTimeCut(iPixel,:);
amp = data.DataAmplitudeCut(iPixel,:);
fitExpSettingsPixel = fitExpSettings;
fitExpSettingsPixel.hParent = hf9;
fitExpSettingsPixel.DoubleExp = true;
fitExpSettingsPixel.SingleExpFitType = 2;
fitExpSettingsPixel.EndOffset = 0;
T = fit_exp(amp*amax, time, fitExpSettingsPixel);


%% -----------------------------------------------------------------------------
% (C) Copyright 2020-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
