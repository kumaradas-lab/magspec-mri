%% Demo "Excitation_Spectra_and_Phase"
% This demonstrates which spectrum and phase is excited by TX pulses of
% different shapes.

%%
% Preparations
LoadSystem
clear TXShape
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 10); % find magnet frequency

% Sequence parameters
Seq.nMeasurement  = 41;                 % number of measurements (e.g. 101)
Seq.nSteadyState  = 4;                  % number of pulses to reach steady state (e.g. 4)
Seq.Bandwidth     = 20e3;               % bandwidth to test in Hz (e.g. 50e3)
Seq.fCenterOffset = 0e3;                % Center frequency offset in Hz (e.g. 0e3)
Seq.tRep          = 200e-3;             % repetition time in seconds (e.g. 300e-3 = 3xT1)
Seq.average       = 1;                  % number of averages (>=1)
Seq.averageBreak  = Seq.tRep;           % time between average measurements (e.g. Seq.tRep)
Seq.useGradSpoil  = 1;                  % select spoiler gradient (0=none, 1=x, 2=y, 3=z)
Seq.GradSpoilTm   = 0.01;               % spoiler gradient strength in T/m
Seq.plotSeqTR     = 0;                  % plot the sequence wrapped at tReps or no sequence plot ([])
Seq.plotSeq       = Seq.useGradSpoil;   % plot the sequence (Seq.useGradSpoil) or no sequence plot ([])


TXShape.usedPulse = 1;                  % choose shape pulse / see switch case below (0, 1, 2, 3)


Seq.tRep          = Seq.tRep*ones(1,floor(Seq.nMeasurement/2)*2+Seq.nSteadyState*2);  % repetition time copy to all tReps
Seq.FrequencyOffsetToTest = [zeros(2,Seq.nSteadyState), ...                     % dummy shots at the center frequency
                            [linspace(Seq.Bandwidth/(Seq.nMeasurement-1),Seq.Bandwidth/2,floor(Seq.nMeasurement/2));... % alternating shots around the center frequency
                            linspace(-Seq.Bandwidth/(Seq.nMeasurement-1),-Seq.Bandwidth/2,floor(Seq.nMeasurement/2))]]; % alternating shots around the center frequency
Seq.FrequencyOffsetToTest = Seq.FrequencyOffsetToTest(:).'+Seq.fCenterOffset;   % put the frequency into a vector
Seq.ntRep = numel(Seq.tRep);


% choose shape for the first pulse
switch TXShape.usedPulse
  case 0 % rect pulse
    TXShape.CenterOffset=   0e-6;                                   % time offset of the center of the segmented pulse
    TXShape.Center=         0e-6;                                   % effective center time of the segmented pulse
    TXShape.Duration=       [HW.tFlip90Def ]+0e-6;                  % durations of pulse segments
    TXShape.Frequency=      [0  ]+HW.fLarmor;                       % frequency of pulse segments
    TXShape.Phase=          [0  ]+0;                                % phase of pulse segments
    TXShape.Amplitude= [HW.TX.AmpDef ]+0*1e-6;                      % amplitude (Tesla) of the B1+ field inside the coil
    TXShape.DurationOfTheSegmentedPulse=sum(TXShape.Duration);      % duration of the whole pulse
    TXShape.Start=cumsum([0;TXShape.Duration(1:end-1)])-TXShape.DurationOfTheSegmentedPulse/2-TXShape.CenterOffset; % calculate the start time

  case 1 % use Pulse_x function
    TXShape.Center=                 0e-6;                           % effective center time of the segmented pulse
    TXShape.BandwidthOfPulse=       10000;                          % bandwidth of the pulse in Hz (e.g. 10000)
    TXShape.FlipAngle=              90;                             % flip angle in degrees (e.g. 90)
    TXShape.maxNumbersOfSegments=   101;                            % maximum number of pulse segments
    TXShape.PulseFunction=          @Pulse_Sinc_2;                  % function of the shaped pulse (e.g. @Pulse_RaisedCos or @Pulse_ press tab)

    %  function [pulseData] = Pulse_RaisedCos(HW, Center, BW, FlipAngleRad, maxpulseCount,  maxLength, Frequency, Phase)
    TXShape=TXShape.PulseFunction(HW,  TXShape.Center,  TXShape.BandwidthOfPulse,  TXShape.FlipAngle/180*pi,  TXShape.maxNumbersOfSegments,  10e3, HW.fLarmor,  0 );

  case 2 % load your own pulse

    load DemoPulseDataAmpPha.mat                                    % load own pulse
    amplitude = amp(:);                                             % amplitude in Hz of B1 fLarmor
    phase = pha(:);                                                 % phase in degrees
    N=size(phase,1);                                                % number of pulse segments
    duration=0.5e-3;                                                % duration of the shaped pulse
    timestep=duration/N;                                            % duration of the pulse segment

    TXShape.CenterOffset=   duration*0.48;                          % time offset of the center of the segmented pulse
    TXShape.Center=         0e-6;                                   % effective center time of the segmented pulse
    TXShape.Duration=       ones(N,1)*timestep;                     % durations of pulse segments
    TXShape.Frequency=      ones(N,1)*HW.fLarmor;                   % frequency of pulse segments
    TXShape.Phase=          phase;                                  % phase of pulse segments
    TXShape.Amplitude=      (amplitude*2*pi/HW.Gamma.H1);           % amplitude (Tesla) of the B1+ field inside the coil
    TXShape.DurationOfTheSegmentedPulse=sum(TXShape.Duration,1);    % duration of the whole pulse
    TXShape.Start=          cumsum([0;TXShape.Duration(1:end-1)],1)-TXShape.DurationOfTheSegmentedPulse/2-TXShape.CenterOffset; % calculate the start time

  case 3 % sinc shaped pulse
    TXShape.NumbersOfZeroCrossings= 6;                              % zero crossings of the sinc function [2,4,6,8,...]
    TXShape.Center=                 0e-6;                           % effective center time of the segmented pulse
    TXShape.BandwidthOfPulse=       10000;                          % bandwidth of the pulse in Hz (e.g. 10000)
    TXShape.FlipAngle=              90;                             % flip angle in degree (e.g. 90)
    TXShape.CenterOffset=           0e-6;                           % time offset of the center of the segmented pulse
    TXShape.maxNumbersOfSegments=   201;                            % maximum number of pulse segments

    TXShape.DurationOfTheSegmentedPulse=TXShape.NumbersOfZeroCrossings/TXShape.BandwidthOfPulse;    % duration of the whole pulse
    if TXShape.maxNumbersOfSegments>TXShape.DurationOfTheSegmentedPulse/(1/HW.TX.fSample*2)         % reduce the number of segments if the pulse is too short
      TXShape.NumbersOfSegments=floor(TXShape.DurationOfTheSegmentedPulse/(1/HW.TX.fSample*2));     % one segment should be at least 1/HW.TX.fSample*2 long (default: 16 ns)
    else
      TXShape.NumbersOfSegments=TXShape.maxNumbersOfSegments;                                       % if the pulse duration is long enough use maxNumbersOfSegments
    end
    TXShape.tShape=      round(linspace(-TXShape.DurationOfTheSegmentedPulse/2+TXShape.DurationOfTheSegmentedPulse/TXShape.NumbersOfSegments,...
      TXShape.DurationOfTheSegmentedPulse/2-TXShape.DurationOfTheSegmentedPulse/TXShape.NumbersOfSegments,...
      TXShape.NumbersOfSegments).'...
      *(HW.TX.fSample/2))/(HW.TX.fSample/2);                                    % generate a timeline rounded to match the TX sampling rate (DA-converter)
    TXShape.Start=       TXShape.tShape-[diff(TXShape.tShape);TXShape.tShape(end)-TXShape.tShape(end-1)]/2-TXShape.CenterOffset;        % calculate the start time so that the tShape time is in the middle of the pulse
    TXShape.Duration=    [diff(TXShape.Start);TXShape.Start(end)-TXShape.Start(end-1)];                       % calculate durations of pulse segments
    TXShape.Frequency=   HW.fLarmor*ones(size(TXShape.Start));                                                                          % frequency of pulse segments
    TXShape.B1Shape=     sinc(TXShape.tShape/TXShape.DurationOfTheSegmentedPulse*TXShape.NumbersOfZeroCrossings);                       % normalized shape of the pulse
    TXShape.Amplitude=   abs(TXShape.B1Shape)*(HW.TX.Amp2FlipPiIn1Sec*(TXShape.FlipAngle/180*pi)/pi)/(TXShape.DurationOfTheSegmentedPulse * mean(TXShape.B1Shape));         % amplitude (Tesla) of the B1+ field in the coil
    TXShape.Phase=       angle(TXShape.B1Shape)/pi*180+0;                              % phase of pulse segments (no negative amplitude is allowed, add 180 to the phase to get an equivalent)


end

TXShape.Start=      TXShape.Start     * ones(1,Seq.ntRep);
TXShape.Duration=   TXShape.Duration  * ones(1,Seq.ntRep);
TXShape.Frequency=  TXShape.Frequency * ones(1,Seq.ntRep) + ones(size(TXShape.Frequency,1),1)*Seq.FrequencyOffsetToTest;
TXShape.Phase=      TXShape.Phase  * ones(1,Seq.ntRep);
TXShape.Amplitude=  TXShape.Amplitude * ones(1,Seq.ntRep);

TXShape=            correct_PulsePhase(TXShape, HW);  % off resonance pulse phase correction


% RF transmission parameters
TX.Start        =   TXShape(1).Start;                             % start times of TX pulse
TX.Duration     =   TXShape(1).Duration;                          % durations of TX pulse
TX.Frequency    =   TXShape(1).Frequency;                         % frequencies of TX pulse
TX.Phase        =   TXShape(1).Phase;                             % phases of TX pulse
TX.Amplitude    =   TXShape(1).Amplitude;                         % amplitudes of TX pulse



% Acquisition parameters
AQ.fSample      =   50e3;                                         % sampling rate of AQ window (e.g. 50e3)
AQ.nSamples     =   max(1e-3*AQ.fSample(1),20);                   % number of samples in AQ window   (e.g. 20)
% acquisition start of AQ window (as close as possible after the first segmented pulse)
AQ.Start        =   TXShape(1).Start(end,:) + TXShape(1).Duration(end,:) + get_DeadTimeTX2RX(HW, AQ.fSample(1));     
AQ.Frequency    =   HW.fLarmor;                                   % frequency of AQ window (e.g. HW.fLarmor)
AQ.Phase        =   0;                                            % phase of AQ window

if Seq.useGradSpoil       % spoiler from end of AQ +1e-3 until Seq.tRep -10e-3 sec
  if Seq.useGradSpoil>3; error('Seq.useGradSpoil set 0 1 2 or 3');end
  Grad(Seq.useGradSpoil).Time    =[   [0;     HW.Grad(1).tRamp] + AQ.Start(1)+AQ.nSamples(1)/AQ.fSample(1)+1e-3;...
    [-HW.Grad(1).tRamp;    0]+Seq.tRep(1)-10e-3]; % time of gradient points in seconds
  Grad(Seq.useGradSpoil).Amp     =[   0;     1;    1;    0;     ]*Seq.GradSpoilTm; % amplitude of gradient points in Tesla/meter
end


% Start measurement
[ Raw, SeqOut, data, data_1D ] = set_sequence(HW, Seq, AQ, TX, Grad);

% zoom into TX pulse
if ~isempty(SeqOut.plotSeqTR)
  hax = get(20, 'Children');
  hax = hax(strcmp(get(hax, 'Type'), 'axes'));
  TXduration = sum(SeqOut.TX(1).Duration(:,1));
  TXcenter = SeqOut.TX(1).Start(1) + TXduration/2;
  set(hax(1), 'XLim', [-TXduration, TXduration] + TXcenter);
end

% Plot results
plot_data_1D(HW, data_1D);

% get Steady State shots
dataSteadyState = data(1).data(:,:,SeqOut.nSteadyState*2:SeqOut.ntRep);  % cut the data of the dummy shots

fMeanAmp = SeqOut.FrequencyOffsetToTest(SeqOut.nSteadyState*2:SeqOut.ntRep).';  % cut the frequency of the dummy shots
[fMeanAmp,sI] = sort(fMeanAmp);                                               % sort frequency and get index

MeanAmp=squeeze(mean(abs(dataSteadyState),1));
MeanPhase=squeeze(mean(unwrap(angle(dataSteadyState)),1));
fOffsetAQ=squeeze(mean(diff(unwrap(angle(dataSteadyState),[],1),[],1),1))/2/pi*SeqOut.AQ(1).fSample(1);
% figure;plot(fOffsetAQ);
MeanfOffset=mean(fOffsetAQ(MeanAmp>max(MeanAmp)/4));
PhaseCorr=(unwrap(angle(squeeze(dataSteadyState).*repmat(exp(-1i*MeanfOffset/(1/2/pi*SeqOut.AQ(1).fSample(1))*(1:SeqOut.AQ(1).nSamples(1))).',1,SeqOut.nMeasurement))));
figure(102)
plot(PhaseCorr);
title PhaseCorr
MeanPhaseCorr=squeeze(mean((angle(squeeze(dataSteadyState).*repmat(exp(1i*MeanfOffset/(1/2/pi*SeqOut.AQ(1).fSample(1))*(1:SeqOut.AQ(1).nSamples(1))).',1,SeqOut.nMeasurement)))));
MeanAmp=MeanAmp(sI);
MeanPhaseCorr=MeanPhaseCorr(sI);

figure(101)
ax(1) = subplot(2,1,1);                           % plot amplitude vs. frequency
plot(ax(1), fMeanAmp, MeanAmp, '-xb', 'LineWidth', 2);
xlabel(ax(1), 'resonance offsets in Hz')
ylabel(ax(1), 'amplitude in Tesla')

ax(2) = subplot(2,1,2);                           % plot phase vs. frequency
plot(ax(2), fMeanAmp, MeanPhaseCorr, '-xr', 'LineWidth', 2);
xlabel(ax(2), 'resonance offsets in Hz')
ylabel(ax(2), 'phase in rad')

linkaxes(ax, 'x');

%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
