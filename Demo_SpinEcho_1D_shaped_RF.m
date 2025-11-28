%% Demo Sequence "Spin Echo 1D with shaped pulses"
% This demo sequence demonstrates how to create a simple spin echo with
% gradients and shaped pulses. The FFT is applied to the encoded spin echo signal which
% shows a profile of the sample.

%% Simple Spin Echo 1D with a shaped RF pulse
% Preparations
LoadSystem;                                 % load system parameters
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 10, [], 1);  % find magnet frequency


% Define parameters required for the measurement
% Parameters used for timing calculations
Seq.tEcho       = 10e-3;                    % echo time ( >=10e-3  )
Seq.p90         = HW.tFlip90Def;            % duration of 1st TX pulse; ( use default value HW.tFlip90Def)
Seq.p180        = HW.tFlip180Def;           % duration of 2nd TX pulse; ( use default value HW.tFlip180Def)

% Sequence parameters
Seq.useGrad     = 1;                        % choose which gradient to be used (1==x, 2==y, 3==z, 0 no gradient)
Seq.SliceGrad   = 70e-3;                    % amplitude of gradient in Tesla/meter during first pulse
Seq.ReadoutGrad = 70e-3;                    % amplitude of gradient in Tesla/meter during echo read-out
Seq.addPhase    = 300;                      % add phase angle too all RF pulses
Seq.shiftGrad   = -0e-6;                    % shift all gradients in time to compensate eddy currents and adjust the position of the spin echo (e.g. -0.1 ms ... +0.1 ms)
Seq.plotSeq     = [Seq.useGrad];            % plot sequence on, plot RF, AQ and Grad(Seq.useGrad) (1==x, 2==y, 3==z, 0 no Grad)
Seq.tRep        = 100e-3;                   % repetition time (>=2*Seq.tEcho and >=T1)
Seq.average     = 10;                       % number of averages (>=1)
Seq.Shape(1).usedPulse = 2;                 % choose shape for the first pulse / see switch case below (0, 1, 2)
Seq.Shape(2).usedPulse = 0;                 % choose shape for the second pulse / see switch case below (0, 1, 2)
Seq.PlotSeqTXB1inHz    = 0;                 % plotSeq TX unit (default 0 Tesla , 1 Hz)


% choose shape for the first Pulse
switch Seq.Shape(1).usedPulse
  case 0 % rect pulse
    Seq.Shape(1).TimeOffsetOfCenterOfPulseSegments=  0e-6;                       % time offset of the center of the segmented pulse
    Seq.Shape(1).DurationsOfPulseSegments=   [Seq.p90 ]+0e-6;                    % durations of pulse segments
    Seq.Shape(1).FrequencyOfPulseSegments=   [0  ]+HW.fLarmor;                   % frequency of pulse segments
    Seq.Shape(1).PhaseOfPulseSegments=       [0  ]+Seq.addPhase;                 % phase of pulse segments
    Seq.Shape(1).B1AmplitudeOfPulseSegments= [HW.TX(1).AmpDef] + 0*1e-6;         % amplitude (Tesla) of the B1+ field inside the coil
    Seq.Shape(1).DurationOfTheSegmentedPulse=sum(Seq.Shape(1).DurationsOfPulseSegments);  % duration of the whole pulse
    Seq.Shape(1).StartOfSegmentedPulses=cumsum([0;Seq.Shape(1).DurationsOfPulseSegments(1:end-1)])-Seq.Shape(1).DurationOfTheSegmentedPulse/2+Seq.Shape(1).TimeOffsetOfCenterOfPulseSegments; % calculate the start time

  case 1 % manually enter the desired values
    Seq.Shape(1).TimeOffsetOfCenterOfPulseSegments=  0e-6;                       % time offset of the center of the segmented pulse
    Seq.Shape(1).DurationsOfPulseSegments=   [10 ;83 ;40 ;20 ;10 ]*1e-6;         % durations of pulse segments
    Seq.Shape(1).FrequencyOfPulseSegments=   [0  ;0  ;0  ;0  ;0  ]+HW.fLarmor;   % frequency of pulse segments
    Seq.Shape(1).PhaseOfPulseSegments=       [0  ;00 ;180;270;0  ]+Seq.addPhase; % phase of pulse segments
    Seq.Shape(1).B1AmplitudeOfPulseSegments= [0 ;3000 ;0; 0;0 ]*2*pi/HW.Gamma.H1;% amplitude (Tesla) of the B1+ field inside the coil, if you use B1 in Hz you have to convert the values using the gyromagnetic ratio of H1. (e.g. 3000 Hz =>  3000*2*pi/HW.Gamma.H1 => 70.46 uT)
    % Seq.Shape(1).B1AmplitudeOfPulseSegments= [10 ;50 ;75; 50;10 ]*1e-6;          % amplitude (Tesla) of the B1+ field in the coil
    Seq.Shape(1).DurationOfTheSegmentedPulse=sum(Seq.Shape(1).DurationsOfPulseSegments);      % duration of the whole pulse
    Seq.Shape(1).StartOfSegmentedPulses=cumsum([0;Seq.Shape(1).DurationsOfPulseSegments(1:end-1)])-Seq.Shape(1).DurationOfTheSegmentedPulse/2+Seq.Shape(1).TimeOffsetOfCenterOfPulseSegments; % calculate the start time

  case 2 % sinc shaped pulse
    Seq.Shape(1).NumbersOfZeroCrossings= 6;                                      % number of zero crossings of the sinc funtion [2,4,6,8,...]
    Seq.Shape(1).BandwidthOfPulse=7000;                                          % bandwidth of the pulse in Hz (e.g. 10000)
    Seq.Shape(1).FlipAngle=50;                                                   % flip angle in degree (e.g. 90)
    Seq.Shape(1).TimeOffsetOfCenterOfPulseSegments=  0e-6;                       % time offset of the center of the segmented pulse
    Seq.Shape(1).maxNumbersOfSegments= 101;                                      % maximum number of pulse segments

    Seq.Shape(1).DurationOfTheSegmentedPulse=Seq.Shape(1).NumbersOfZeroCrossings/Seq.Shape(1).BandwidthOfPulse; % duration of the whole pulse
    if Seq.Shape(1).maxNumbersOfSegments>Seq.Shape(1).DurationOfTheSegmentedPulse/(1/HW.TX(1).fSample*2)        % reduce the number of segments if the pulse is too short
      Seq.Shape(1).NumbersOfSegments=floor(Seq.Shape(1).DurationOfTheSegmentedPulse/(1/HW.TX(1).fSample*2));  % one segment should be at least 1/HW.TX.fSample*2 long (default: 16 ns)
    else
      Seq.Shape(1).NumbersOfSegments = Seq.Shape(1).maxNumbersOfSegments;        % if the pulse duration is long enough use maxNumbersOfSegments
    end
    Seq.Shape(1).tShape=round(linspace(-Seq.Shape(1).DurationOfTheSegmentedPulse/2+Seq.Shape(1).DurationOfTheSegmentedPulse/Seq.Shape(1).NumbersOfSegments,...
                                        Seq.Shape(1).DurationOfTheSegmentedPulse/2-Seq.Shape(1).DurationOfTheSegmentedPulse/Seq.Shape(1).NumbersOfSegments,...
                                        Seq.Shape(1).NumbersOfSegments).'...
                                        *(HW.TX(1).fSample/2))/(HW.TX(1).fSample/2);  % generate a timeline rounded to match the TX sampling rate (DA-converter)
    Seq.Shape(1).StartOfSegmentedPulses=Seq.Shape(1).tShape-[diff(Seq.Shape(1).tShape);Seq.Shape(1).tShape(end)-Seq.Shape(1).tShape(end-1)]/2+Seq.Shape(1).TimeOffsetOfCenterOfPulseSegments;       % calculate the start time so that the tShapte time is in the middle of the pulse
    Seq.Shape(1).DurationsOfPulseSegments=   [diff(Seq.Shape(1).StartOfSegmentedPulses);Seq.Shape(1).StartOfSegmentedPulses(end)-Seq.Shape(1).StartOfSegmentedPulses(end-1)];                       % calculate durations of pulse segments
    Seq.Shape(1).FrequencyOfPulseSegments=   zeros(Seq.Shape(1).NumbersOfSegments,1)+HW.fLarmor;                                                                                                    % frequency of pulse segments
    Seq.Shape(1).B1Shape=                    sinc(Seq.Shape(1).tShape/Seq.Shape(1).DurationOfTheSegmentedPulse*Seq.Shape(1).NumbersOfZeroCrossings);                                                % normalized shape of the pulse
    Seq.Shape(1).B1AmplitudeOfPulseSegments=   abs(Seq.Shape(1).B1Shape)*(HW.TX(1).Amp2FlipPiIn1Sec*(Seq.Shape(1).FlipAngle/180*pi)/pi)/sum(Seq.Shape(1).DurationsOfPulseSegments .* Seq.Shape(1).B1Shape,1);  % amplitude (Tesla) of the B1+ field in the coil
    Seq.Shape(1).PhaseOfPulseSegments=       angle(Seq.Shape(1).B1Shape)/pi*180+0+Seq.addPhase;  % phase of pulse segments (no negative amplitude is allowed, so you have to add 180 degrees to the phase to get an equivalent)

end

% choose shape for the second pulse
switch Seq.Shape(2).usedPulse
  case 0 % rect pulse
    Seq.Shape(2).TimeOffsetOfCenterOfPulseSegments=  Seq.tEcho/2 + 0e-6;         % time offset of the center of the segmented pulse
    Seq.Shape(2).DurationsOfPulseSegments=   [Seq.p180 ]+0e-6;                   % durations of pulse segments
    Seq.Shape(2).FrequencyOfPulseSegments=   [0  ]+HW.fLarmor;                   % frequency of pulse segments
    Seq.Shape(2).PhaseOfPulseSegments=       [90  ]+Seq.addPhase;                % phase of pulse segments
    Seq.Shape(2).B1AmplitudeOfPulseSegments= [HW.TX(1).AmpDef ]+0*1e-6;          % amplitude (Tesla) of the B1+ field in the coil
    Seq.Shape(2).DurationOfTheSegmentedPulse=sum(Seq.Shape(2).DurationsOfPulseSegments);  % duration of the whole pulse
    Seq.Shape(2).StartOfSegmentedPulses=cumsum([0;Seq.Shape(2).DurationsOfPulseSegments(1:end-1)])-Seq.Shape(2).DurationOfTheSegmentedPulse/2+Seq.Shape(2).TimeOffsetOfCenterOfPulseSegments;  % calculate the Start Time

  case 1 % put in the values by hand
    Seq.Shape(2).TimeOffsetOfCenterOfPulseSegments=  Seq.tEcho/2 + 0e-6;         % time offset of the center of the segmented pulse
    Seq.Shape(2).DurationsOfPulseSegments=   [10 ;20 ;40 ;20 ;10 ]*1e-6;         % durations of pulse segments
    Seq.Shape(2).FrequencyOfPulseSegments=   [0  ;0  ;0  ;0  ;0  ]+HW.fLarmor;   % frequency of pulse segments
    Seq.Shape(2).PhaseOfPulseSegments=       [0  ;90 ;180;270;0  ]+Seq.addPhase; % phase of pulse segments
    Seq.Shape(2).B1AmplitudeOfPulseSegments= [10 ;50 ;75; 50;10 ]*1e-6;          % amplitude (Tesla) of the B1+ field in the coil
    Seq.Shape(2).DurationOfTheSegmentedPulse=sum(Seq.Shape(2).DurationsOfPulseSegments);  % duration of the whole pulse
    Seq.Shape(2).StartOfSegmentedPulses=cumsum([0;Seq.Shape(2).DurationsOfPulseSegments(1:end-1)])-Seq.Shape(2).DurationOfTheSegmentedPulse/2+Seq.Shape(2).TimeOffsetOfCenterOfPulseSegments; % calculate the Start Time

  case 2 % sinc shaped pulse
    Seq.Shape(2).NumbersOfZeroCrossings= 6;                                      % number ot times the sinc function crosses the zero [2,4,6,8,...]
    Seq.Shape(2).BandwidthOfPulse=5000;                                          % bandwidth of the pulse in Hz (e.g. 5000)
    Seq.Shape(2).FlipAngle=180;                                                  % flip angle in degree (e.g. 180)
    Seq.Shape(2).TimeOffsetOfCenterOfPulseSegments=  Seq.tEcho/2 + 0e-6;         % time offset of the center of the segmented pulse
    Seq.Shape(2).maxNumbersOfSegments= 101;                                      % maximum numbers of pulse segments

    Seq.Shape(2).DurationOfTheSegmentedPulse=Seq.Shape(2).NumbersOfZeroCrossings/Seq.Shape(2).BandwidthOfPulse; % duration of the whole pulse
    if Seq.Shape(2).maxNumbersOfSegments>Seq.Shape(2).DurationOfTheSegmentedPulse/(1/HW.TX(1).fSample*2)        % reduce the number of segments if the pulse is too short
      Seq.Shape(2).NumbersOfSegments=floor(Seq.Shape(2).DurationOfTheSegmentedPulse/(1/HW.TX(1).fSample*2));  % one segment should be at least 1/HW.TX.fSample*2 long (default: 16 ns)
    else
      Seq.Shape(2).NumbersOfSegments=Seq.Shape(2).maxNumbersOfSegments;         % if the pulse is long enough use maxNumbersOfSegments
    end
    Seq.Shape(2).tShape=round(linspace(-Seq.Shape(2).DurationOfTheSegmentedPulse/2+Seq.Shape(2).DurationOfTheSegmentedPulse/Seq.Shape(2).NumbersOfSegments,...
                                        Seq.Shape(2).DurationOfTheSegmentedPulse/2-Seq.Shape(2).DurationOfTheSegmentedPulse/Seq.Shape(2).NumbersOfSegments,...
                                        Seq.Shape(2).NumbersOfSegments).'...
                                        *(HW.TX(1).fSample/2))/(HW.TX(1).fSample/2);  % generate a timeline rounded to match the TX sampling rate (DA-converter)
    Seq.Shape(2).StartOfSegmentedPulses=Seq.Shape(2).tShape-[diff(Seq.Shape(2).tShape);Seq.Shape(2).tShape(end)-Seq.Shape(2).tShape(end-1)]/2+Seq.Shape(2).TimeOffsetOfCenterOfPulseSegments;   % calculate the start time so that the tShapte time is in the middle of the pulse
    Seq.Shape(2).DurationsOfPulseSegments= [diff(Seq.Shape(2).StartOfSegmentedPulses);Seq.Shape(2).StartOfSegmentedPulses(end)-Seq.Shape(2).StartOfSegmentedPulses(end-1)];                     % calculate durations of pulse segments
    Seq.Shape(2).FrequencyOfPulseSegments= zeros(Seq.Shape(2).NumbersOfSegments,1)+HW.fLarmor;                                                                                                  % frequency of pulse segments
    Seq.Shape(2).B1Shape=                    sinc(Seq.Shape(2).tShape/Seq.Shape(2).DurationOfTheSegmentedPulse*Seq.Shape(2).NumbersOfZeroCrossings);                                            % normalized shape of the pulse
    Seq.Shape(2).B1AmplitudeOfPulseSegments= abs(Seq.Shape(2).B1Shape)*(HW.TX(1).Amp2FlipPiIn1Sec*(Seq.Shape(2).FlipAngle/180*pi)/pi)/sum(Seq.Shape(2).DurationsOfPulseSegments .* Seq.Shape(2).B1Shape,1);  % amplitude (Tesla) of the B1+ field in the coil
    Seq.Shape(2).PhaseOfPulseSegments=       angle(Seq.Shape(2).B1Shape)/pi*180+90+Seq.addPhase;                            % phase of pulse segments (no negative smplitude is allowed, so you have to add 180 degrees to the phase to get an equivalent)

end

if Seq.PlotSeqTXB1inHz  % only used in plotSeq and plotSeqTR for plotting
  HW.TX(1).AmplitudeUnitScale = (2*pi/HW.Gamma.H1);
  HW.TX(1).AmplitudeUnit = 'Hz';
  HW.TX(1).AmplitudeName = 'TX B1';
end

% RF transmission parameters
TX.Start      = [ Seq.Shape(1).StartOfSegmentedPulses; ...        % start times of 1st TX pulses
                  Seq.Shape(2).StartOfSegmentedPulses];           % start times of 2nd TX pulses
TX.Duration   = [ Seq.Shape(1).DurationsOfPulseSegments; ...      % durations of 1st TX pulses
                  Seq.Shape(2).DurationsOfPulseSegments];         % durations of 2nd TX pulses
TX.Frequency  = [ Seq.Shape(1).FrequencyOfPulseSegments; ...      % frequencies of 1st TX pulses
                  Seq.Shape(2).FrequencyOfPulseSegments];         % frequencies of 2nd TX pulses
TX.Phase      = [ Seq.Shape(1).PhaseOfPulseSegments; ...          % phases of 1st TX pulses
                  Seq.Shape(2).PhaseOfPulseSegments];             % phases of 2nd TX pulses
TX.Amplitude  = [ Seq.Shape(1).B1AmplitudeOfPulseSegments; ...    % amplitudes of 1st TX pulses
                  Seq.Shape(2).B1AmplitudeOfPulseSegments];       % amplitudes of 2nd TX pulses


% Acquisition parameters
AQ.fSample    = [ 50e3; ...                   % sampling rate of 1st AQ window
                  30e3];                      % sampling rate of 2st AQ window
AQ.nSamples   = [ 128; ...                    % number of samples in 1st AQ window
                  128];                       % number of samples in 2nd AQ window
AQ.Start      = [ TX.Start(length(Seq.Shape(1).DurationsOfPulseSegments))+TX.Duration(length(Seq.Shape(1).DurationsOfPulseSegments))+get_DeadTimeTX2RX(HW,AQ.fSample(1)); ... % acquisition start of 1st AQ window (as close as possible after the first segmented pulse)
                  Seq.tEcho-(AQ.nSamples(2)+1)/AQ.fSample(2)/2]; % acquisition start of 2nd AQ window (the center of the window is at Seq.tEcho)
AQ.Frequency  = [ HW.fLarmor; ...             % frequency of 1st AQ window
                  HW.fLarmor];                % frequency of 2nd AQ window
AQ.Phase      = [ 0; ...                      % phase of 1st AQ window
                  145];                       % phase of 2nd AQ window

% Gradient along x, y or z (1==x, 2==y, 3==z, 0 no Grad)
if Seq.useGrad
  Grad(Seq.useGrad).Time  = [ [-0.6; -0.5;  0.5; 0.6 ]*1e-3;...
                              [0.7 ;  0.9;  3.2; 3.4 ]*1e-3;...
                              [-2.6; -2.4;  2.4; 2.6 ]*1e-3+Seq.tEcho]; % time of gradient points in seconds
  Grad(Seq.useGrad).Amp   = [ [0;    1;    1;    0]*Seq.SliceGrad;...
                              [0;    1;    1;    0]*(Seq.ReadoutGrad-Seq.SliceGrad*1.1/2.5/2);...
                              [0;    1;    1;    0]*Seq.ReadoutGrad  ]; % amplitude of gradient points in Tesla/meter
  if (AQ.Start(2)<Grad(Seq.useGrad).Time(7))
    error('AQ.Start(2) to early, decrease AQ.nSamples(2), or increase AQ.fSample(2)');
  end
  HW.Grad(1).TimeDelay = HW.Grad(1).TimeDelay - Seq.shiftGrad;
else
  for t = 1:4
    Grad(t).Time = NaN;  % gradient not used
    Grad(t).Amp = 0;     % no amplitude
  end
end

% Start measurement
[Raw, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);


% Plot results
plot_data_1D(HW, data_1D);

% Plot FFT
hf = figure(1);
hax = subplot(3,1,1, 'Parent', hf);
plot(hax, ...
  (data(1).time_of_tRep(1:AQ.nSamples(2),2))*1000, abs(data(1).data(1:AQ.nSamples(2),2))*1e9, ...
  (data(1).time_of_tRep(1:AQ.nSamples(2),2))*1000, real(data(1).data(1:AQ.nSamples(2),2))*1e9, ...
  (data(1).time_of_tRep(1:AQ.nSamples(2),2))*1000, imag(data(1).data(1:AQ.nSamples(2),2))*1e9);
xlabel(hax, 'time in ms');
ylabel(hax, 'amplitude in nT');
grid(hax, 'on');
legend(hax, 'abs', 'real', 'imag');

hax = subplot(3,1,2, 'Parent', hf);
plot(hax, ...
  (data(1).f_fft1_data(1:AQ.nSamples(2),2)-HW.fLarmor)/1000, abs(data(1).fft1_data(1:AQ.nSamples(2),2)), ...
  (data(1).f_fft1_data(1:AQ.nSamples(2),2)-HW.fLarmor)/1000, real(data(1).fft1_data(1:AQ.nSamples(2),2)), ...
  (data(1).f_fft1_data(1:AQ.nSamples(2),2)-HW.fLarmor)/1000, imag(data(1).fft1_data(1:AQ.nSamples(2),2)));
xlabel(hax, 'offset frequency in kHz');
ylabel(hax, 'amplitude');
grid(hax, 'on');

hax = subplot(3,1,3, 'Parent', hf);
plot(hax, ...
  (data(1).f_fft1_data(1:AQ.nSamples(2),2)-HW.fLarmor)/HW.Gamma.H1*2*pi/Seq.ReadoutGrad*1000, abs(data(1).fft1_data(1:AQ.nSamples(2),2)),...
  (data(1).f_fft1_data(1:AQ.nSamples(2),2)-HW.fLarmor)/HW.Gamma.H1*2*pi/Seq.ReadoutGrad*1000, real(data(1).fft1_data(1:AQ.nSamples(2),2)),...
  (data(1).f_fft1_data(1:AQ.nSamples(2),2)-HW.fLarmor)/HW.Gamma.H1*2*pi/Seq.ReadoutGrad*1000, imag(data(1).fft1_data(1:AQ.nSamples(2),2)))
xlabel(hax, 'distance in mm')
ylabel(hax, 'amplitude')
grid(hax, 'on');

%% -----------------------------------------------------------------------------
% (C) Copyright 2011-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
