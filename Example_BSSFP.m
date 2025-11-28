%% Balanced Steady-State Free Precession
%
% This sequence can be used e.g. to record a MRI "video".
% Do proper shimming with Demo_Auto_ParameterSearch for best results!
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

LoadSystem;                                               % load system parameters (reset to default: HW Seq AQ TX Grad)

Seq.T1Estimated = 500e-3;
Seq.AQSlice(1).nImages = 40;                              % number of acquired images
Seq.AQSlice(1).TurboBreak = max(2, Seq.T1Estimated*3);    % pause between two loop averages in seconds ([]= fast as possible)

Seq.LoopSaveAllData = 1;                                  % save all data for each image (tEcho)
Seq.LoopPlotAverages = 0;                                 % do not plot average of tEchoes
Seq.AQSlice(1).plotkSpace = 0;

Seq.Find_Frequency_interval = 1;                          % run Find_Frequency every 1000 seconds or after every loop (if measurement duration > 1000 sec)

Seq.Loops = 1;                                            % number of loop averages 1...

Seq.tEcho = 3.8e-3;                                       % Echo time in seconds e.g. 5e-3


% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).nRead = 40;                                % number of pixels in read, if nRead>1 nPhase(1)=1
Seq.AQSlice(1).nPhase(1) = 1;                             % number of pixels in phase(1) (3D)
Seq.AQSlice(1).nPhase(2) = 40;                            % number of pixels in phase(2) (2D)
Seq.AQSlice(1).nPhase(3) = 1;                             % number of pixels in phase(3) (CSI)
Seq.AQSlice(1).HzPerPixMin = 1000;                        % bandwidth per pixel in Hz (1/HzPerPixMin= duration of AQ)
Seq.AQSlice(1).sizeRead = 10e-3;                          % size in read direction in meter
Seq.AQSlice(1).sizePhase(2) = 10e-3;                      % size in phase(2) direction in meter
Seq.AQSlice(1).thickness = 4e-3;                          % slice thickness of excitation pulse in meter
Seq.AQSlice(1).thicknessInversion = Seq.AQSlice(1).thickness;  % slice thickness of inversion pulse in meter (bandwidth!)

% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice(1).PhaseOS(1) = 1;                            % oversampling phase(1)  1...
Seq.AQSlice(1).PhaseOS(2) = 2;                            % oversampling phase(2)  1...
Seq.AQSlice(1).PhaseOS(3) = 1;                            % oversampling phase(3)  1...
Seq.AQSlice(1).PhaseOS(Seq.AQSlice(1).nPhase==1) = 1;     % reset PhaseOS(x) to 1 if nPhase(x) == 1
Seq.AQSlice(1).ReadOS = 1;

% % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nPhase = prod(Seq.AQSlice(1).nPhase .* Seq.AQSlice(1).PhaseOS);
Seq.AQSlice(1).oddEvenEchoes = 0;                         % acquire separate images consisting of only odd or even Echoes respectively
Seq.AQSlice(1).phaseCycling = 0;                          % acquire separate images with the inversion pulses rotated by 180 degrees (supress Echoes of inversion pulses)
Seq.AQSlice(1).TurboBlocks = 1;

% % Orientation in space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'xz');  % Top-down image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'zx');  % Bottom-up image (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yx');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))
% Seq.AQSlice = get_AlphaPhiTheta(Seq.AQSlice, 'yz');  % image encoding directions (slice/phase(1), phase(2), read/phase(3))

% % Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Seq.plotSeq = 1:3;                                        % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
Seq.LoopPlot = 1;                                         % plot every loop
Seq.plot3D = 0;
Seq.AQSlice(1).plotkSpace = 1;                            % plot k-space
Seq.AQSlice(1).plotImage = 1;                             % plot image
Seq.AQSlice(1).plotPhase = 1;                             % plot phase of k-space and/or image
Seq.AQSlice(1).ZeroFillWindowSize = 1.4;                  % zero fill window size (high k-space values are damped by a cos^2 law)
Seq.AQSlice(1).ZeroFillFactor = 2;                        % zero fill resolution factor
Seq.AQSlice(1).partsAverage = 1;

Seq.MaxGradAmpSlice = 0.1;                                % limit slice gradient strength
Seq.CorrectSliceRephase = 0*double(Seq.AQSlice(1).thickness<=0.015); % correct SliceGradTimeIntegralRephaseOffset
Seq.CorrectPhaseRephase = 0;                              % correct PhaseGradTimeIntegralRephaseOffset
Seq.CorrectReadRephase = 0;                               % correct ReadGradTimeIntegralOffset
Seq.AQSlice(1).SpoilFactor = [0, 0.00, 0.00];             % spoiler resolution factor (slice or phase(1), phasae 2, read or phase 3)
Seq.AQSlice(1).SpoilLengthFactor = 1;                     % spoil Length Factor


% BSSFP
Seq.AQSlice(1).excitationFlipAngle = 45;
Seq.AQSlice(1).inversionFlipAngle = 90;
Seq.AQSlice(1).readOutPhaseIncrement = 180;
Seq.AQSlice(1).inversionPhaseIncrement = 180;
Seq.AQSlice(1).excitationPulse = @Pulse_Sinc_3_Hamming;   % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
Seq.AQSlice(1).inversionPulse = @Pulse_Sinc_3_Hamming;    % inversion pulse function
Seq.AQSlice(1).readoutPhaseInversion = true;              % read out with phase of inversion pulses
% Seq.AQSlice(1).excitationPhase = 90;


Seq.SteadyState_PreShots90 = 0;
Seq.SteadyState_PreShots180 = 8;

[SeqLoop, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave);  % run measurement


if 0
  %% Set new Zero Fill parameter
  SeqLoop.AQSlice(1).ZeroFillWindowSize = 1.11;               % zero fill window size (high k-space values are damped by a cos^2 law)
  SeqLoop.AQSlice(1).ZeroFillFactor = 2;                      % zero fill resolution factor
  for t = SeqLoop.nEchos:-1:1
    SeqLoop.AQSlice(t).ZeroFillWindowSize = SeqLoop.AQSlice(1).ZeroFillWindowSize;
    SeqLoop.AQSlice(t).ZeroFillFactor = SeqLoop.AQSlice(1).ZeroFillFactor;
    % get data of echo
    SeqLoop.data = get_kSpaceAndImage(SeqLoop.data,SeqLoop.AQSlice(t));
    if SeqLoop.nEchos > 1
      if ~isfield(SeqLoop.data, 'ImageZAll')
        SeqLoop.data.ImageZAll = nan(size(SeqLoop.data.ImageZ,1), size(SeqLoop.data.ImageZ,2), size(SeqLoop.data.ImageZ,3), size(SeqLoop.data.ImageZ,4), SeqLoop.nEchos);
      end
      SeqLoop.data.ImageZAll(:,:,:,:,t) = SeqLoop.data.ImageZ;
    end
    % data.Image(max(abs(data.Image(:)))/4>abs(data.Image(:))) = nan;
  end
end


%% Collect and process data from SeqLoop
clear data
% get images
data.ImageZAll = SeqLoop.data.ImageZ(:,:,:,:,:,1);
% phase correction with first image
data.ImageZAll = bsxfun(@times, data.ImageZAll, exp(-1i*(0+angle(data.ImageZAll(:,:,:,:,1)))));
% get k-space of all images
data.kSpaceOsRawAll = SeqLoop.data.kSpaceOsRaw(:,:,:,:,:,1);
% phase correction with first k-space
data.kSpaceOsRawAll = bsxfun(@times, data.kSpaceOsRawAll, exp(-1i*(0+angle(data.kSpaceOsRawAll(:,:,:,:,1)))));

amax = max(abs(data.ImageZAll(:)));

% permute dimension with tEchos to front in amplitude
data.DataAmplitude = permute(real(data.ImageZAll), [5,1,2,3,4])/amax;
data.DataTime = SeqLoop.data.tImageZ;

% sliceomatic(1, squeeze(abs(data.ImageZAll))) % Plot T2 image stack
% Plot some images of the T2 image stack
t = 0;
hf = figure(2);
clf(2,'reset')
minmax = [-max(abs(data.ImageZAll(:))),max(abs(data.ImageZAll(:)))]/amax;
amaxk = max(abs(data.kSpaceOsRawAll(:)));
minmaxk = [-max(abs(data.kSpaceOsRawAll(:))), max(abs(data.kSpaceOsRawAll(:)))]/amaxk;

numAll = unique(round(linspace(1, SeqLoop.AQSlice(1).nImages, 6)));

if Seq.AQSlice(1).plotkSpace, k=2; else k=1; end
for num = numAll
    t = t+1;
    minmax2 = [-max(max(max(max(max(abs(data.ImageZAll(:,:,:,:,num)))))))-eps(amax), ...
      max(max(max(max(max(abs(data.ImageZAll(:,:,:,:,num)))))))+eps(amax)]/amax;
    % image amplitude
    hax = subplot(numel(numAll),5*k,1+(t-1)*5*k);
    imagesc(squeeze(abs(data.ImageZAll(:,:,:,:,num))).'/amax, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmax2);
    title(hax, ['abs ' num2str(1/minmax2(2),'%.1f')]);
    ylabel(hax, ['TE' num2str(num) '=' num2str(data.DataTime(num)),' s']);
    % image real
    hax = subplot(numel(numAll),5*k,2+(t-1)*5*k, 'Parent', hf);
    imagesc(squeeze(real(data.ImageZAll(:,:,:,:,num))).'/amax, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmax2);
    title(hax, ['real ' num2str(1/minmax2(2),'%.1f')]);
    % image imag
    hax = subplot(numel(numAll),5*k,3+(t-1)*5*k, 'Parent', hf);
    imagesc(squeeze(imag(data.ImageZAll(:,:,:,:,num))).'/amax, 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', minmax2);
    title(hax, ['imag ' num2str(1/minmax2(2),'%.1f')]);
    % image phase
    hax = subplot(numel(numAll),5*k,4+(t-1)*5*k, 'Parent', hf);
    imagesc(squeeze(angle(data.ImageZAll(:,:,:,:,num))).', 'Parent', hax);
    set(hax, 'YDir', 'normal', 'CLim', [-pi,pi]);
    title(hax, 'phase');
    % image amplitude with color limits of first row
    hax = subplot(numel(numAll),5*k,5+(t-1)*5*k, 'Parent', hf);
    imagesc(squeeze(abs(data.ImageZAll(:,:,:,:,num))).'/amax, 'Parent', hax);
    set(hax,'YDir','normal', 'CLim', minmax);
    title(hax, ['abs ' num2str(1/minmax2(2),'%.1f') ', const. CLim']);
    if k==2
      minmaxk2 = [-max(max(max(max(max(abs(data.kSpaceOsRawAll(:,:,:,:,num)))))))-eps(amaxk), ...
        max(max(max(max(max(abs(data.kSpaceOsRawAll(:,:,:,:,num)))))))+eps(amaxk)]/amaxk;
      % kSpace amplitude
      hax = subplot(numel(numAll),5*k,1+(t-1)*5*k+5, 'Parent', hf);
      imagesc(log10(squeeze(abs(data.kSpaceOsRawAll(:,:,:,:,num))).'/amaxk), 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', [-3 log10(minmaxk2(2))]);
      title(hax, ['abs ' num2str(1/minmaxk2(2),'%.1f')]);
      % kSpace real
      hax = subplot(numel(numAll),5*k,2+(t-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(real(data.kSpaceOsRawAll(:,:,:,:,num))).'/amaxk, 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', minmaxk2);
      title(hax, ['real ' num2str(1/minmaxk2(2),'%.1f')]);

      hax = subplot(numel(numAll),5*k,3+(t-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(imag(data.kSpaceOsRawAll(:,:,:,:,num))).'/amaxk, 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', minmaxk2);
      title(hax, ['imag ' num2str(1/minmaxk2(2),'%.1f')]);

      hax = subplot(numel(numAll),5*k,4+(t-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(angle(data.kSpaceOsRawAll(:,:,:,:,num))).', 'Parent', hax);
      set(hax, 'YDir', 'normal', 'CLim', [-pi,pi]);
      title(hax, 'phase');

      hax = subplot(numel(numAll),5*k,5+(t-1)*5*k+5, 'Parent', hf);
      imagesc(squeeze(abs(data.kSpaceOsRawAll(:,:,:,:,num))).'/amaxk, 'Parent', hax);
      set(hax,'YDir','normal', 'CLim', minmaxk);
      title(hax, 'abs 1.00');
    end

end
