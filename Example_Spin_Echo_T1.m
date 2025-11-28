%% T1 map
%
% Saturation or inversion pulse with increasing delay to a following
% Turbo Spin Echo imaging block.

LoadSystem;                                                 % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.Recovery = 'Inversion';                                 % 'Saturation' or 'Inversion'
Seq.T1Estimated = 3000e-3;                                  % estimated T1 of sample in seconds, e.g. 3 for water
Seq.T2Estimated = 2000e-3;                                  % estimated T1 of sample in seconds, e.g. 3 for water

% pause between two loops in seconds ([]= fast as possible)
if strcmp(Seq.Recovery, 'Inversion')
  Seq.LoopsBreak = max(Seq.T1Estimated*3, 0.8);
else
  Seq.LoopsBreak = 0.8;
end

HW.FindFrequencyPause = max(Seq.LoopsBreak-1, 1);           % pause after freequency sweep in seconds

Seq.Function_Prepare_Measurement = @prepare_Recovery;
Seq.tFlipStart = Seq.T1Estimated/5;                         % flip wait time at start, if zero use minimum  e.g. 10e-3          (T1/10)
Seq.tFlipEnd = Seq.T1Estimated*5;                           % flip wait time at end,                        e.g. 500e-3         (5*T1)
Seq.tFlipLog = 1;                                           % increase wait time logarithmically (1) or linearly (0)
Seq.LoopSaveAllData = 1;
Seq.LoopPlotAverages = 0;
Seq.LoopPlotLastAverage = 0;
Seq.LoopSeqPlot = 1;

Seq.Loops = 7;                                              % number of loops, i.e., different preparation times

Seq.Find_Frequency_interval = 0;                            % time between two frequency searches in seconds
% Seq.plotSeqTR = [1:3];                                      % plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
Seq.plotSeq = 1:3;                                          % plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
% Seq.plotSeqStart = 10;                                      % plot sequence start with tRep x
% Seq.plotSeqEnd = 11;                                        % plot sequence stop with tRep x
Seq.plot_iLaplaceMean = 0;                                  % plot mean of iLaplace T1 spectrum per pixel
Seq.plot_iLaplace3D = 0;                                    % plot iLaplace T1 fit

Seq.tEcho = 4e-3;                                           % Echo time in seconds e.g. 5e-3

% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 24;                                  % number of pixels in read, if nRead>1 nPhase(1)=1
Seq.AQSlice(1).nPhase(1) = 1;                               % number of pixels in phase(1)
Seq.AQSlice(1).nPhase(2) = 24;                              % number of pixels in phase(2)
Seq.AQSlice(1).HzPerPixMin = 1000;                          % bandwith per pixel in Hz (1/HzPerPixMin= duration of AQ)
Seq.AQSlice(1).sizeRead = 0.0096;                           % image size in read in meter (for CSI set to 1e12)
Seq.AQSlice(1).sizePhase(1) = 0.0096;                       % image size in phase(1) in meter
Seq.AQSlice(1).sizePhase(2) = 0.0096;                       % image size in phase(2) in meter
dim_image = sum([Seq.AQSlice(1).nPhase, Seq.AQSlice(1).nRead] > 1);
if dim_image < 3
  % select slice only for 2d image
  Seq.AQSlice(1).thickness = 0.005;                         % slice thickness in meter, used for 2D and 3D! ([] for no slice)
end
Seq.AQSlice(1).excitationPulse = @Pulse_RaisedCos;          % excitation pulse function (type "Pulse_" and press tab for selection of pulses)
Seq.AQSlice(1).inversionPulse = @Pulse_Rect_Composite180;   % inversion pulse function

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.AQSlice(1).ReadOS = 16;                                 % oversampling read ([] for automatic, recommended >=2) 1...
Seq.AQSlice(1).PhaseOS(1) = 1;                              % oversampling phase(1)  1...
Seq.AQSlice(1).PhaseOS(2) = 2;                              % oversampling phase(2)  1...
Seq.AQSlice(1).PhaseOS(3) = 1;                              % oversampling phase(3)  1... (set to 1 if nPhase(3)=1;
Seq.AQSlice(1).PhaseOS(Seq.AQSlice(1).nPhase==1) = 1;       % reset PhaseOS(x) to 1 if nPhase(x) == 1

% % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.AQSlice(1).TurboFactor = 1;                             % number of image k-lines per excitation
% Seq.AQSlice(1).TurboFactor = prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS)/1;  % number of image k-lines per excitation
% Seq.tRep = max(Seq.LoopsBreak,(Seq.tFlipEnd+50e-3)*double(Seq.AQSlice(1).TurboFactor==1));  % repetition time in seconds only used if Seq.TurboBreak=[], auto Seq.tRep=[] and Seq.AQSlice(1).TurboBreak= x sec;
% Seq.nPreLoops = double(Seq.AQSlice(1).TurboFactor>1);
% Seq.AQSlice(1).TurboBreak = Seq.tRep;                       % break between last echo and next excitation

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.LoopPlot = 0;                                         % plot every loop
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
Seq.AQSlice(1).Center2OriginImage = [0.00,0.000,0.000];     % vector from center of the image to the origin in image coordinate system

% % Some corrections    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.AQSlice(1).SpoilLengthFactor = 1;                       % spoil length factor
% Seq.AQSlice(1).sizePhaseSpoil = [1,1,1]*0.0005;%1e12;       % spoil size
if dim_image < 3
  % larger spoiler in slice direction for 2d images
  Seq.AQSlice(1).SpoilFactor = [10, 1, 1];
else
  % equivalent spoiler amplitudes in all image directions for 3d images
  Seq.AQSlice(1).SpoilFactor = [1, 1, 1];
end
Seq.MaxGradAmpSlice = 0.05;                                 % limit slice gradient strength
% Seq.SteadyState_PreShots180 = 0;

%% Run experiment
[SeqLoop, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);

if 0
  %% Re-calculate k-space and image (with differing settings)
  SeqLoop.AQSlice(1).ZeroFillWindowSize = 1.;               % zero fill window size (high k-space values are damped by a cos^2 law)
  SeqLoop.AQSlice(1).ZeroFillFactor = 2;                    % zero fill resolution factor
  SeqLoop.dataLoop(1).StartSequenceTime = SeqLoop.data.StartSequenceTime;
  SeqLoop.dataLoop(1).fCenter = SeqLoop.data.fCenter;
  SeqLoop.dataLoop(1).kTicks = SeqLoop.data.kTicks;
  SeqLoop.dataLoop(1).Ticks = SeqLoop.data.Ticks;
  SeqLoop.dataLoop(1).RoICutOff = SeqLoop.data.RoICutOff;
  SeqLoop.dataLoop(1).RoI = SeqLoop.data.RoI;
  for t = 1:SeqLoop.Loops
    [SeqLoop.dataLoop(t)] = get_kSpaceAndImage(SeqLoop.dataLoop(t), SeqLoop.AQSlice(1));
    [SeqLoop.dataLoop(t)] = plot_kSpaceAndImage(SeqLoop.dataLoop(t), SeqLoop.AQSlice(1));
  end
  [SeqLoop.data] = get_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
  [SeqLoop.data] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
end


%% evaluation
clear data
data.ImageZAll = cat(5, SeqLoop.dataLoop(SeqLoop.nPreLoops+(1:SeqLoop.nFlips)).ImageZ);
% data.ImageZAll = data.ImageZAll.*exp(-1i*(0+angle((squeeze(data.ImageZAll(Seq.AQSlice(1).nRead*Seq.AQSlice.ZeroFillFactor/2,1,Seq.AQSlice(1).nPhase(2)*Seq.AQSlice.ZeroFillFactor/2,1,1))))));
% data.ImageZAll = bsxfun(@times, data.ImageZAll,exp(-1i*(0+angle(data.ImageZAll(:,:,:,:,1)))));
data.ImageZAll = bsxfun(@times, data.ImageZAll, exp(-1i*(0+angle(data.ImageZAll(:,:,:,:,end)))));

dim_image = sum([SeqLoop.AQSlice(1).nPhase, SeqLoop.AQSlice(1).nRead] > 1);
if dim_image < 3
  % display a selection of 2d images for some preparation times
  iAxes = 0;
  hf = figure(2);
  amax = max(abs(data.ImageZAll(:)));
  minmax = [-max(abs(data.ImageZAll(:))), max(abs(data.ImageZAll(:)))]/amax;
  for num = max(round([1,(SeqLoop.nFlips-1)/5+1,(SeqLoop.nFlips-1)*2/5+1,(SeqLoop.nFlips-1)*3/5+1,(SeqLoop.nFlips-1)*4/5+1,SeqLoop.nFlips]), 1)
    iAxes = iAxes+1;
    minmax2 = [-max(max(max(max(max(abs(data.ImageZAll(:,:,:,:,num))))))), ...
               max(max(max(max(max(abs(data.ImageZAll(:,:,:,:,num)))))))]/amax;
    hax = subplot(6,5,1+(iAxes-1)*5, 'Parent', hf);
    imagesc(squeeze(abs(data.ImageZAll(:,:,:,:,num))).'/amax, 'Parent', hax, minmax2);
    set(hax, 'YDir', 'normal');
    title(hax, 'abs');
    ylabel(hax, ['t_e ' num2str(SeqLoop.tFlip(SeqLoop.nPreLoops+num)) ' s']);

    hax = subplot(6,5,2+(iAxes-1)*5, 'Parent', hf);
    imagesc(squeeze(real(data.ImageZAll(:,:,:,:,num))).'/amax, 'Parent', hax, minmax2)
    set(hax,'YDir','normal')
    title(hax, 'real')

    hax = subplot(6,5,3+(iAxes-1)*5, 'Parent', hf);
    imagesc(squeeze(imag(data.ImageZAll(:,:,:,:,num))).'/amax, 'Parent', hax, minmax2)
    set(hax,'YDir','normal')
    title(hax, 'imag')

    hax = subplot(6,5,4+(iAxes-1)*5, 'Parent', hf);
    imagesc(squeeze(angle(data.ImageZAll(:,:,:,:,num))).', 'Parent', hax, [-pi,pi])
    set(hax,'YDir','normal')
    title(hax, 'angle')

    hax = subplot(6,5,5+(iAxes-1)*5, 'Parent', hf);
    imagesc(squeeze(real(data.ImageZAll(:,:,:,:,num))).'/amax, 'Parent', hax, minmax)
    set(hax,'YDir','normal')
    title(hax, 'real norm')
  end
end

data.DataAmplitude = permute(real(data.ImageZAll),[5,1,2,3,4]);
data.DataTime = SeqLoop.tFlip(SeqLoop.nPreLoops+(1:SeqLoop.nFlips));


%% Single-exponential fit of amplitude in each pixel
nT1Estimated = 6;

n = min(size(data.DataAmplitude,1), ceil(nT1Estimated*SeqLoop.T1Estimated/diff(data.DataTime(1:2))));
if dim_image < 3
  % 2d image
  RoI = SeqLoop.data.RoI(:,:,1);
  SeqLoop.data.RoIAll = bsxfun(@times, reshape(RoI, 1, size(RoI,1), 1, size(RoI,2)), [ones(n,1);nan(size(data.DataAmplitude,1)-n,1)]);
else
  % 3d image
  RoI = SeqLoop.data.RoI(:,:,:,1);
  SeqLoop.data.RoIAll = bsxfun(@times, reshape(RoI, [1, size(RoI)]), [ones(n,1);nan(size(data.DataAmplitude,1)-n,1)]);
end
data.DataTimeCut = repmat(data.DataTime(1:n), sum(~isnan(RoI(:))), 1);
data.DataAmplitudeCut = reshape(data.DataAmplitude(~isnan(SeqLoop.data.RoIAll)), n, []).';

T1map = nan(size(RoI));
T1 = nan(sum(~isnan(RoI(:))), 1);
clear fitExpSettings
fitExpSettings.EndOffset = true;
fitExpSettings.DoubleExp = false;
fitExpSettings.SingleExpFitType = 2;
fitExpSettings.RingFilter = 0;
% t=round(sum(~isnan(SeqLoop.data.RoI(:))/2))+28;
for iPixel = 1:sum(~isnan(RoI(:)))
  x = data.DataTimeCut(iPixel,:);
  y = data.DataAmplitudeCut(iPixel,:);
  [T, fitExpSettings] = fit_exp(y, x, fitExpSettings);
  T1(iPixel) = T.tau;
end
T1map(~isnan(RoI(:))) = T1;

hf10 = figure(10);
if dim_image < 3
  % 2d image
  if ~isempty(getappdata(hf10, 'sliceomatic'))
    % previous 3d results are still open
    clf(hf10);
  end
  ax10(1) = subplot(1,1,1, 'Parent', hf10);
  imagesc(SeqLoop.data.Ticks(1).ReadZ, SeqLoop.data.Ticks(2).PhaseZ, T1map.', 'Parent', ax10(1));
  set(ax10(1), 'CLim', [0, SeqLoop.T1Estimated*2]);
  set(ax10(1), 'YDir', 'normal');
  title(ax10(1), ['T1 in s, mean = ' num2str(mean(T1)) ' s, STD = ' num2str(std(T1)) ' s']);
  xlabel(ax10(1), 'Read in m');
  ylabel(ax10(1), 'Phase in m');
  colorbar('peer', ax10(1));
  set(ax10(1), 'DataAspectRatio', [1 1 1]);
else
  % 3d image
  hsl = sliceomatic(hf10, T1map, SeqLoop.data.Ticks(1).ReadZ, SeqLoop.data.Ticks(1).PhaseZ, SeqLoop.data.Ticks(2).PhaseZ);
  set(hsl, 'CLim', [0, SeqLoop.T1Estimated*2]);
  title(hsl.hAxes, ['T1 in s, mean = ' num2str(mean(T1)) ' s, STD = ' num2str(std(T1)) ' s']);
  xlabel(hsl.hAxes, [SeqLoop.AQSlice(1).ReadCartesianAxis{1}, 'Read in m']);
  ylabel(hsl.hAxes, [SeqLoop.AQSlice(1).PhaseCartesianAxis{1}, 'Phase(1) in m']);
  zlabel(hsl.hAxes, [SeqLoop.AQSlice(1).PhaseCartesianAxis{2}, 'Phase(2) in m']);
  title(hsl.GetSliderX(), SeqLoop.AQSlice(1).ReadCartesianAxis{1});
  title(hsl.GetSliderY(), SeqLoop.AQSlice(1).PhaseCartesianAxis{1});
  title(hsl.GetSliderZ(), SeqLoop.AQSlice(1).PhaseCartesianAxis{2});
  if isa(hsl, 'sliceomatic') && ...
      isempty([hsl.GetAllSlicesPosX(); hsl.GetAllSlicesPosY(); hsl.GetAllSlicesPosZ(); hsl.GetAllIsoValues()])
    % Note: The axes labels don't necessarily correspond to the axes "orientation".
    hsl.AddSliceX(0);
  end
end
set(hf10, 'Name', 'T1 (Single-Exponential Fit)');

%% Plot data at singular pixel
hf9 = figure(9);
set(hf9, 'Name', 'Single pixel data')

iPixel = round(sum(~isnan(RoI(:))/2))+0;
x = data.DataTimeCut(iPixel,:);
y = data.DataAmplitudeCut(iPixel,:);
fitExpSettingsPixel = fitExpSettings;
fitExpSettingsPixel.hParent = hf9;
T = fit_exp(y, x, fitExpSettingsPixel);


%% Calculate inverse Laplace transform
if SeqLoop.plot_iLaplaceMean || (SeqLoop.plot_iLaplace3D && dim_image < 3)
  % SeqOut.iLaplace1D.SpectrumTimeStart = [];
  % SeqOut.iLaplace1D.SpectrumTimeEnd = [];
  % SeqOut.iLaplace1D.SpectrumTimeStartCut = [];
  % SeqOut.iLaplace1D.SpectrumTimeEndCut = [];
  SeqLoop.iLaplace1D.FullScaleAmplitude = [];
  SeqLoop.iLaplace1D.Problem = SeqLoop.Recovery; % 'Inversion' or 'Saturation'
  SeqLoop.iLaplace1D.SpectrumTimeStart    = data.DataTime(1)/2;
  SeqLoop.iLaplace1D.SpectrumTimeEnd      = data.DataTime(end)*2;
  SeqLoop.iLaplace1D.SpectrumTimeStartCut = data.DataTime(1)/1;
  SeqLoop.iLaplace1D.SpectrumTimeEndCut   = data.DataTime(end)*1;
  SeqLoop.iLaplace1D.nSpectrum = 100;
  SeqLoop.iLaplace1D.Plot  =0;
  SeqLoop.iLaplace1D.QualityFactor = 100;
  [data, SeqLoop] = get_iLaplace1D(data, SeqLoop);
  % data.iLaplace1D.SpectrumAmplitude(data.iLaplace1D.SpectrumAmplitude>1)=0;
  % data.iLaplace1D.SpectrumAmplitude(data.iLaplace1D.SpectrumAmplitude<0)=0;
end


%% Plot mean T1 values of inverse Laplace transform
if SeqLoop.plot_iLaplaceMean
  data.iLaplace1D.SpectrumTimeMean = bsxfun(@times, ...
    sum(bsxfun(@times, data.iLaplace1D.SpectrumTime, data.iLaplace1D.SpectrumAmplitude(:,:,:,:,:)),1), ...
    1/sum(data.iLaplace1D.SpectrumAmplitude,1));
  % data.iLaplace1D.SpectrumTimeMean=bsxfun(@times,sum(bsxfun(@times,data.iLaplace1D.SpectrumTime(2:end-1),data.iLaplace1D.SpectrumAmplitude(2:end-1,:,:,:,:)),1),1/sum(data.iLaplace1D.SpectrumAmplitude,1));
  hf4 = figure(4);
  if dim_image < 3
    % 2d image
    if ~isempty(getappdata(hf4, 'sliceomatic'))
      % previous 3d results are still open
      clf(hf4);
    end
    ax4(1) = subplot(1,1,1, 'Parent', hf4);
    imagesc(SeqLoop.data.Ticks(1).ReadZ,SeqLoop.data.Ticks(2).PhaseZ, (squeeze(data.iLaplace1D.SpectrumTimeMean) .* SeqLoop.data.RoI).', 'Parent', ax4(1))
    set(ax4(1), 'CLim', [0,SeqLoop.T1Estimated*2])
    set(ax4(1), 'YDir', 'normal')
    title(ax4(1), 'mean T1 in s')
    xlabel(ax4(1), 'Read in m')
    ylabel(ax4(1), 'Phase in m')
    colorbar('peer', ax4(1));
    set(ax4(1), 'DataAspectRatio', [1 1 1])
  else
    % 3d image
    hsl = sliceomatic(hf4, permute(squeeze(data.iLaplace1D.SpectrumTimeMean), SeqLoop.data.PermuteOrder) .* SeqLoop.data.RoI, ...
      SeqLoop.data.Ticks(1).ReadZ, SeqLoop.data.Ticks(1).PhaseZ, SeqLoop.data.Ticks(2).PhaseZ);
    set(hsl, 'CLim', [0, SeqLoop.T1Estimated*2]);
    title(hsl.hAxes, 'mean T1 in s')
    xlabel(hsl.hAxes, [SeqLoop.AQSlice(1).ReadCartesianAxis{1}, 'Read in m']);
    ylabel(hsl.hAxes, [SeqLoop.AQSlice(1).PhaseCartesianAxis{1}, 'Phase(1) in m']);
    zlabel(hsl.hAxes, [SeqLoop.AQSlice(1).PhaseCartesianAxis{2}, 'Phase(2) in m']);
    title(hsl.GetSliderX(), SeqLoop.AQSlice(1).ReadCartesianAxis{1});
    title(hsl.GetSliderY(), SeqLoop.AQSlice(1).PhaseCartesianAxis{1});
    title(hsl.GetSliderZ(), SeqLoop.AQSlice(1).PhaseCartesianAxis{2});
    if isa(hsl, 'sliceomatic') && ...
        isempty([hsl.GetAllSlicesPosX(); hsl.GetAllSlicesPosY(); hsl.GetAllSlicesPosZ(); hsl.GetAllIsoValues()])
      % Note: The axes labels don't necessarily correspond to the axes "orientation".
      hsl.AddSliceX(0);
    end
  end
  set(hf4, 'Name', 'Mean T1 of Inverse Laplace Transform')
end


%% Plot iLaplace T1 fit
if SeqLoop.plot_iLaplace3D && dim_image < 3
  hf = figure(3);
  close(hf);
  hf = figure(3);
  hsl = sliceomatic(hf, cat(1, squeeze(bsxfun(@times, SeqLoop.iLaplace1D.FullScaleAmplitude,data.iLaplace1D.SpectrumAmplitude)), ...
                               sum(squeeze(bsxfun(@times, SeqLoop.iLaplace1D.FullScaleAmplitude, data.iLaplace1D.SpectrumAmplitude)), 1) ...
                                  / max(max(sum(squeeze(bsxfun(@times, SeqLoop.iLaplace1D.FullScaleAmplitude, data.iLaplace1D.SpectrumAmplitude)), 1))) ...
                                  * max(max(max(squeeze(bsxfun(@times, SeqLoop.iLaplace1D.FullScaleAmplitude, data.iLaplace1D.SpectrumAmplitude))))) ...
                            ), ...
                      SeqLoop.data.Ticks(1).ReadZ, ...
                      [data.iLaplace1D.SpectrumTime;data.iLaplace1D.SpectrumTime(end)+diff(data.iLaplace1D.SpectrumTime(end-1:end))].', ...
                      SeqLoop.data.Ticks(2).PhaseZ);
  set(hsl, 'YScale', 'log');
  set(hsl, 'DataAspectRatioMode', 'auto');
  set(hsl, 'PlotBoxAspectRatio', [1,1,1]);
  % set(hsl, 'PlotBoxAspectRatioMode', 'manual');
end


%% -----------------------------------------------------------------------------
% (C) Copyright 2016-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
