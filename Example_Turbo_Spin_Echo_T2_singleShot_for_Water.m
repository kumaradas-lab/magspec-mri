%% Example script to measure a 2d T2 map with a single shot Turbo Spin Echo
%
% Do proper shimming with Demo_Auto_ParameterSearch for best results.

LoadSystem;                                               % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.T2Estimated = 2.150;
Seq.tEchoTrain  = Seq.T2Estimated*10;                     % duration of echo Train
Seq.LoopsBreak  = max(2, Seq.T2Estimated*5);              % pause between two loop averages in seconds ([]= fast as possible)

Seq.LoopSaveAllData = 1;                                  % save all data for each image (tEcho)
Seq.LoopPlotAverages = 0;                                 % do not plot average of tEchoes

HW.Grad(1).tRamp = 0.2e-3;                                % increase raise and fall time of gradients
HW.Grad(1).tEC   = 0.1e-3;                                % adjust time delay of gradient pulse

HW.FindFrequencyPause = Seq.T2Estimated*5;                % increase pause after frequency search
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 1);       % find magnet frequency
Seq.Find_Frequency_interval = 1000;                       % Run Find_Frequency every 1000 seconds or after every loop (if measurement duration > 1000 sec)

Seq.Loops = 1;                                            % number of loop averages 1...

Seq.tEcho = 12e-3;                                        % echo time in seconds e.g. 5e-3

Seq.SteadyState_PreShots180=0;
% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 20;                                % number of pixels in read, if nRead>1 nPhase(1)=1
Seq.AQSlice(1).nPhase(1) = 1;                             % number of pixels in phase(1) (3D)
Seq.AQSlice(1).nPhase(2) = 20*1;                          % number of pixels in phase(2) (2D)
Seq.AQSlice(1).nPhase(3) = 1;                             % number of pixels in phase(3) (CSI)
Seq.AQSlice(1).HzPerPixMin = 1000;                        % bandwidth per pixel in Hz (1/HzPerPixMin= duration of AQ)
Seq.AQSlice(1).sizeRead = 10e-3;                          % size in read direction in meter
Seq.AQSlice(1).sizePhase(2) = 10e-3;                      % size in phase(2) direction in meter
Seq.AQSlice(1).thickness = 0.005;                         % slice thickness in meter % slice not working properly yet
Seq.AQSlice(1).thicknessInversion = Seq.AQSlice(1).thickness*1.5; % slice thickness of inversion pulses in meter
Seq.AQSlice(1).excitationPulse = @Pulse_Sinc_3_Hamming;   % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.AQSlice(1).inversionPulse = @Pulse_Sinc_3_Hamming;    % inversion pulse function

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(1) = 1;                            % oversampling phase(1)  1...
Seq.AQSlice(1).PhaseOS(2) = 2;                            % oversampling phase(2)  1...
Seq.AQSlice(1).PhaseOS(3) = 1;                            % oversampling phase(3)  1...
Seq.AQSlice(1).PhaseOS(Seq.AQSlice(1).nPhase==1) = 1;     % reset PhaseOS(x) to 1 if nPhase(x) == 1

% % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).oddEvenEchoes = 1;                         % acquire separate images consisting of only odd or even Echoes respectively
Seq.AQSlice(1).phaseCycling = 0;                          % acquire separate images with the inversion pulses rotated by 180 degrees
Seq.AQSlice(1).nImages = ceil(Seq.tEchoTrain/(Seq.tEcho*prod(Seq.AQSlice(1).nPhase.*Seq.AQSlice(1).PhaseOS))/(Seq.AQSlice(1).oddEvenEchoes+1));  % number of acquired images
Seq.AQSlice(1).TurboFactor = prod(Seq.AQSlice(1).nPhase.*Seq.AQSlice(1).PhaseOS)/1;      % number of image k-lines per excitation
Seq.AQSlice(1).TurboBreak = Seq.LoopsBreak;               % break between last echo and next excitation
Seq.tRep = Seq.LoopsBreak;                                % Pause between two loop averages in seconds ([]= fast as possible)

% % Spoiling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).SpoilFactor = [0, 0.00, 0.00];             % spoiler resolution factor (slice or phase 1, phase 2, read or phase 3)
Seq.AQSlice(1).SpoilFactorInversionBlock = 8;

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xz');  % Top-down image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zx');  % Bottom-up image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeq = 1:3;                                        % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                         % plot every loop
Seq.plot3D = 0;
Seq.AQSlice(1).plotkSpace = 0;                            % plot k-space
Seq.AQSlice(1).plotImage = 1;                             % plot image
Seq.AQSlice(1).plotPhase = 1;                             % plot phase of k-space or image
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;                  % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 2;                        % zero fill resolution factor

Seq.MaxGradAmpSlice = 0.02;                               % limit slice gradient strength
Seq.CorrectSliceRephase = 0*double(Seq.AQSlice(1).thickness<=0.015); % correct SliceGradTimeIntegralRephaseOffset
Seq.CorrectPhaseRephase = 0;                              % correct PhaseGradTimeIntegralRephaseOffset
Seq.CorrectReadRephase = 0;                               % correct ReadGradTimeIntegralOffset
Seq.AQSlice(1).SpoilFactor = [0, 0.00, 0.00];             % spoiler resolution factor (slice or phase 1, phase 2, read or phase 3)
Seq.AQSlice(1).SpoilLengthFactor = 1;                     % spoil length factor


[SeqLoop, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);  % run measurement


if 0
  %% Set new Zero Fill parameter
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
t = 0;
hf = figure(2);
clf(2,'reset')
minmax = [-max(abs(data.ImageZ(:))),max(abs(data.ImageZ(:)))]/amax;
amaxk = max(abs(data.kSpaceOsRaw(:)));
minmaxk = [-max(abs(data.kSpaceOsRaw(:))), max(abs(data.kSpaceOsRaw(:)))]/amaxk;
numAll = ceil([1, SeqLoop.AQSlice(1).nImages/81, SeqLoop.AQSlice(1).nImages/27,SeqLoop.AQSlice(1).nImages/9,SeqLoop.AQSlice(1).nImages/3,SeqLoop.AQSlice(1).nImages-1]);
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
if Seq.AQSlice(1).plotkSpace, k=2; else k=1; end
for num = numAll
  t = t+1;
  minmax2 = [-max(max(max(max(max(abs(data.ImageZ(:,:,:,:,num)))))))-eps(amax), ...
    max(max(max(max(max(abs(data.ImageZ(:,:,:,:,num)))))))+eps(amax)]/amax;
  % image amplitude
  hax = subplot(numel(numAll),5*k,1+(t-1)*5*k);
  imagesc(squeeze(abs(data.ImageZ(:,:,:,:,num))).'/amax, 'Parent', hax);
  set(hax, 'YDir', 'normal', 'CLim', minmax2);
  title(hax, ['abs ' num2str(1/minmax2(2),'%.1f')]);
  ylabel(hax, ['TE' num2str(num) '=' num2str(data.DataTime(num)),' s']);
  % image real
  hax = subplot(numel(numAll),5*k,2+(t-1)*5*k, 'Parent', hf);
  imagesc(squeeze(real(data.ImageZ(:,:,:,:,num))).'/amax, 'Parent', hax);
  set(hax, 'YDir', 'normal', 'CLim', minmax2);
  title(hax, ['real ' num2str(1/minmax2(2),'%.1f')]);
  % image imag
  hax = subplot(numel(numAll),5*k,3+(t-1)*5*k, 'Parent', hf);
  imagesc(squeeze(imag(data.ImageZ(:,:,:,:,num))).'/amax, 'Parent', hax);
  set(hax, 'YDir', 'normal', 'CLim', minmax2);
  title(hax, ['imag ' num2str(1/minmax2(2),'%.1f')]);
  % image phase
  hax = subplot(numel(numAll),5*k,4+(t-1)*5*k, 'Parent', hf);
  imagesc(squeeze(angle(data.ImageZ(:,:,:,:,num))).', 'Parent', hax);
  set(hax, 'YDir', 'normal', 'CLim', [-pi,pi]);
  title(hax, 'phase');
  % image amplitude with color limits of first row
  hax = subplot(numel(numAll),5*k,5+(t-1)*5*k, 'Parent', hf);
  imagesc(squeeze(abs(data.ImageZ(:,:,:,:,num))).'/amax, 'Parent', hax);
  set(hax,'YDir','normal', 'CLim', minmax);
  title(hax, ['abs ' num2str(1/minmax2(2),'%.1f') ', const. CLim']);
  if k==2
    minmaxk2 = [-max(max(max(max(max(abs(data.kSpaceOsRaw(:,:,:,:,num)))))))-eps(amaxk), ...
      max(max(max(max(max(abs(data.kSpaceOsRaw(:,:,:,:,num)))))))+eps(amaxk)]/amaxk;
    % kSpace amplitude
    hax = subplot(numel(numAll),5*k,1+(t-1)*5*k+5, 'Parent', hf);
    imagesc(squeeze(abs(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmaxk2);
    title(hax, ['abs ' num2str(1/minmaxk2(2),'%.1f')]);
    % kSpace real
    hax = subplot(numel(numAll),5*k,2+(t-1)*5*k+5, 'Parent', hf);
    imagesc(squeeze(real(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmaxk2);
    title(hax, ['real ' num2str(1/minmaxk2(2),'%.1f')]);

    hax = subplot(numel(numAll),5*k,3+(t-1)*5*k+5, 'Parent', hf);
    imagesc(squeeze(imag(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmaxk2);
    title(hax, ['imag ' num2str(1/minmaxk2(2),'%.1f')]);

    hax = subplot(numel(numAll),5*k,4+(t-1)*5*k+5, 'Parent', hf);
    imagesc(squeeze(angle(data.kSpaceOsRaw(:,:,:,:,num))).', 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', [-pi,pi]);
    title(hax, 'phase');

    hax = subplot(numel(numAll),5*k,5+(t-1)*5*k+5, 'Parent', hf);
    imagesc(squeeze(abs(data.kSpaceOsRaw(:,:,:,:,num))).'/amaxk, 'Parent', hax);
    set(hax,'YDir','normal', 'CLim', minmaxk);
    title(hax, 'abs 1.00');
  end

end
drawnow

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

%%

% SeqLoop.iLaplace1D.Problem = SeqLoop.Recovery; % 'Inversion' or 'Saturation'
SeqLoop.iLaplace1D.SpectrumTimeStart = data.DataTime(1)/2;
SeqLoop.iLaplace1D.SpectrumTimeEnd  = data.DataTime(end)*2;
SeqLoop.iLaplace1D.SpectrumTimeStartCut = SeqLoop.iLaplace1D.SpectrumTimeStart*2;%;
SeqLoop.iLaplace1D.SpectrumTimeEndCut = SeqLoop.iLaplace1D.SpectrumTimeEnd*0.99;%;
SeqLoop.iLaplace1D.nSpectrum = 200;
SeqLoop.iLaplace1D.Plot = 0;
SeqLoop.iLaplace1D.QualityFactor = 20;

[data, SeqLoop] = get_iLaplace1D(data, SeqLoop);

% Plot mean T2 values
data.iLaplace1D.SpectrumTimeMean = bsxfun(@times, sum(bsxfun(@times, data.iLaplace1D.SpectrumTime, data.iLaplace1D.SpectrumAmplitude(:,:,:,:,:)), 1), 1/sum(data.iLaplace1D.SpectrumAmplitude, 1));
hf4 = figure(4);
ax4(1) = subplot(1,1,1, 'Parent', hf4);
RoI = SeqLoop.data.RoI(:,:,1);
imagesc(SeqLoop.data.Ticks(1).ReadZ, SeqLoop.data.Ticks(2).PhaseZ,squeeze(data.iLaplace1D.SpectrumTimeMean).' .* RoI.', 'Parent', ax4(1));
set(ax4(1), 'CLim', [0,SeqLoop.T2Estimated*2]);
set(ax4(1), 'YDir', 'normal');
title(ax4(1), ...
  ['T2 / s, mean = ' num2str(mean(data.iLaplace1D.SpectrumTimeMean(~isnan(RoI)), 'omitnan')) ', ', ...
  'STD = ' num2str(std(data.iLaplace1D.SpectrumTimeMean(~isnan(RoI)), 'omitnan'))])
xlabel(ax4(1), 'Read in m');
ylabel(ax4(1), 'Phase in m');
colorbar('peer', ax4(1));
set(ax4(1), 'DataAspectRatio', [1 1 1]);
set(hf4, 'Name', 'T2 (Mean of Inverse Laplace Transformation)', 'NumberTitle', 'off');

% Plot iLaplace T2 fit
if SeqLoop.plot3D
  hf = figure(3);
  close(hf);
  hf = figure(3);
  hsl = sliceomatic(hf,cat(1, squeeze(bsxfun(@times, SeqLoop.iLaplace1D.FullScaleAmplitude,data.iLaplace1D.SpectrumAmplitude)),...
                              sum(squeeze(bsxfun(@times, SeqLoop.iLaplace1D.FullScaleAmplitude,data.iLaplace1D.SpectrumAmplitude)),1)...
                                  /max(max(sum(squeeze(bsxfun(@times, SeqLoop.iLaplace1D.FullScaleAmplitude,data.iLaplace1D.SpectrumAmplitude)),1)))...
                                  *max(max(max(squeeze(bsxfun(@times, SeqLoop.iLaplace1D.FullScaleAmplitude,data.iLaplace1D.SpectrumAmplitude)))))...
                          ),...
                      SeqLoop.data.Ticks(1).ReadZ,...
                      [data.iLaplace1D.SpectrumTime;data.iLaplace1D.SpectrumTime(end)+diff(data.iLaplace1D.SpectrumTime(end-1:end))].',...
                      SeqLoop.data.Ticks(2).PhaseZ...
                  );
  set(hsl, 'YScale', 'log')
  set(hsl, 'DataAspectRatioMode', 'auto')
  set(hsl, 'PlotBoxAspectRatio', [1,1,1])
end

%% Single-exponential fit of amplitude in each pixel
nT2Estimated = 10;

n = min(size(data.DataAmplitude,1), ceil(nT2Estimated*SeqLoop.T2Estimated/diff(data.DataTime(1:2))));
RoI = SeqLoop.data.RoI(:,:,1);
SeqLoop.data.RoIAll = bsxfun(@times, reshape(RoI, 1, size(RoI,1), 1, size(RoI,2)), [ones(n,1);nan(size(data.DataAmplitude,1)-n,1)]);

data.DataTimeCut = repmat(data.DataTime(1:n), 1, sum(~isnan(RoI(:))));
data.DataAmplitudeCut = abs(reshape(data.DataAmplitude(~isnan(SeqLoop.data.RoIAll)), n, []));


T2imag = nan(size(RoI));
T2 = nan(sum(~isnan(RoI(:))), 1);
for t = 1:sum(~isnan(RoI(:)))
  x = data.DataTimeCut(:,t);
  y = data.DataAmplitudeCut(:,t);
  w = exp(((y.*exp(-x/(SeqLoop.T2Estimated))).^0.5)/exp(1))-1; % good
  w = w/w(1)*y(1); % normalize weights for plot
  W = diag(w);

  yLog = log(y); % linearize data
  X = [ones(length(x),1) x];
  % X = [y x];
  % b = (X\yLog)
  b = ((W*X)\(w.*yLog)); % solve linearized equation
  % yfit = exp((b(1)+x.*b(2)));

  T2(t) = -1/b(2);
end
T2imag(~isnan(RoI(:))) = T2;

hf10 = figure(10);
ax10(1) = subplot(1,1,1, 'Parent', hf10);
imagesc(SeqLoop.data.Ticks(1).ReadZ, SeqLoop.data.Ticks(2).PhaseZ, T2imag.', 'Parent', ax10(1));
set(ax10(1), 'CLim', [0,SeqLoop.T2Estimated*2]);
set(ax10(1), 'YDir', 'normal');
title(ax10(1), ['T2 / s, mean = ' num2str(mean(T2)) ', STD = ' num2str(std(T2))]);
xlabel(ax10(1), 'Read in m');
ylabel(ax10(1), 'Phase in m');
colorbar('peer',ax10(1));
set(ax10(1), 'DataAspectRatio', [1 1 1]);
set(hf10, 'Name', 'T2 (Single-Exponential Fit)', 'NumberTitle', 'off')


%% Plot data at singular pixel
hf9 = figure(9);
set(hf9, 'Name', 'Singular pixel data')

t = round(sum(~isnan(RoI(:))/2))+-30;
x = data.DataTimeCut(:,t);
y = data.DataAmplitudeCut(:,t);
w = exp(((y.*exp(-x/(SeqLoop.T2Estimated))).^0.5)/exp(1))-1; % good
w = w/w(1)*y(1); % normalize weights for plot
W = diag(w);

yLog = log(y); % linearize data
X = [ones(length(x),1) x];
% X = [y x];
% b = (X\yLog)
b = ((W*X)\(w.*yLog)); % solve linearized equation
yfit = exp((b(1)+x.*b(2)));

T2s = -1/b(2);

hax = subplot(2,1,1, 'Parent', hf9);
semilogy(x,y, x,yfit, x,w, 'Parent', hax);
xlabel(hax, 'time in s');
ylabel(hax, 'amplitude (norm)');
grid(hax, 'on');
hax = subplot(2,1,2, 'Parent', hf9);
plot(x,y, x,yfit, x,w, 'Parent', hax);
xlabel(hax, 'time in s')
ylabel(hax, 'amplitude (norm)')
legend(hax, {'measured data', 'linearized fit', 'weights for lin. fit'})

% sum((y-yfit).^2)^0.5

grid(hax, 'on');


%% -----------------------------------------------------------------------------
% (C) Copyright 2019-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
