%% Demo Sequence "Auto Parameter Search"
% automatically find the correct frequency, RF 90 and 180 degree pulses and
% shim values; values will be appended to "User\MagnetShimCal.m" and
% "PaUout2AmplitudeCal.m"
% Use sample which is only in the homogeneous region of the rf coil.
% T1 of about 100 ms for fast and good results.


%% Auto Parameter Search
% Preparations
LoadSystem;                           % load system parameters

T1 = 0.15;  % approximate T1 of sample in seconds
HW.FindFrequencyPause = T1*2;  % pause duration after frequency search in seconds

%% Find_Frequency_Sweep( HW, mySave, maxtime, span, doplot, tPulse90)
% parameters:
%  HW, mySave; these parameters are created automatically by LoadSystem
LoadSystem;                           % load system parameters
HW.FindFrequencyPause = T1*2;

minTime               = 1;            % minimum time in seconds since the last Find_Frequency; i.e. if set to 10, find freq will only be run, if 10 seconds have passed since the last search
span                  = 200e3;        % searching span of frequency around fLarmor (HW.B0*HW.FindFrequencyGamma/2/pi)
doplot                = 1;            % plot sequence and data
tPulse90              = HW.tFlip90Def/3;  % duration of the 90 degrees pulse ([] uses HW.FindFrequencySweep.tPulse90 or HW.tFlip90)
nMeasurements         = 21;           % number of measurements (TX and AQ trains)
nSamples              = 100;          % number of samples in each measurement
HW.FindFrequencySweep.fSample = 20e3;
HW.FindFrequencySweep.fOffsetFIDsStdMaxValue = 10;
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, minTime, span, doplot, tPulse90, nMeasurements, nSamples);


%% Find_Shim(HW, mySave, maxtime, doplot, iterations, tEcho, T1, ShimStart, ShimStep)
% parameters:
%  HW, mySave; these parameters are created automatically by LoadSystem
LoadSystem;                           % load system parameters
HW.FindFrequencyPause = T1*2;
% HW.FindFrequencySweep.fOffsetFIDsStdMaxValue = 200
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);

minTime               = 10;           % minimum time in seconds since the last Find_Frequency_Sweep.
doplot                = 1;            % plot sequence and data
Seq.iterations        = 100;          % number of iterations used for shim search
Seq.tEcho             = min(T1*4,1);  % echo time
Seq.ShimStart         = [HW.MagnetShim];  % start with these 1x4 shim values ([] for defaults)
Seq.ShimStep          = 1e-4;         % 1x4 vector or scalar with step widths ([] for defaults)
Seq.nPreLoop          = 5;            % number of loops before fminsearch (T1)
Seq.nEchos            = 0;            % use FID or Echo for shimming (0 for Fid 1,2,3... for nEcho)
Seq.use_nEchos_Frequency = 1;         % find frequency of nEchos
Seq.RepetitionTime    = max([1/3, Seq.tEcho*(Seq.nEchos+0.5), T1/10, (Seq.nEchos>0)*T1]);  % repetition time: ~estimated T1 value of the used sample * 3 e.g. 0.5 s (Seq.RepetitionTime=3*T1)
Seq.excitationFlipAngle = acosd(exp(-Seq.RepetitionTime/T1));

clear SliceSelect
Seq.useSliceSelect    = 0;            % use slice gradient
SliceSelect.thickness = 0.008;        % slice thickness
SliceSelect.MaxGradAmpSlice = 0.005;  % maximum gradient amplitude during slice excitation
% SliceSelect.alfa      = 0.0*pi;       % rotation around x axis
% SliceSelect.phi       = 0.0*pi;       % rotation around y axis
% SliceSelect.theta     = 0.5*pi;       % rotation around z axis
SliceSelect = get_AlphaPhiTheta(SliceSelect, 'yzx');  % imaging encoding directions (slice/phase(1), phase(2), read/phase(3))
if Seq.useSliceSelect
  Seq.plotSeq = 1:3;
end


[HW, mySave] = Find_Shim(HW, mySave, minTime, doplot, Seq, SliceSelect);

LoadSystem;                           % load new parameters from file.
HW.FindFrequencyPause=T1*2;

disp(['New shim values: ' num2str(HW.MagnetShim) ' T/m.']);

% search frequency with new shim values
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);


%% Find excitation pulse
% It applies an excitation pulse and acquires the FID after some perpendicular pulses.
LoadSystem;                           % load system parameters
HW.FindFrequencyPause = T1*2;
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);

Seq.T1 = T1;
Seq.tExcitation = 100e-6;  % fix duration of 90 degrees pulse in seconds
Seq.setExcitation = 'Time';  % 'Time' 'Amp' 'PaUout' or []
Seq.FixedExcitationPulseLength = true;

% Seq.PaUoutExcitation = 1;
% Seq.setExcitation = 'PaUout';  % 'Time', 'Amp', 'PaUout', or []
% Seq.FixedExcitationPulseLength = false;

Seq.SearchFlipDeg = 90*1;
Seq.AQPhaseOffsetDegSaveLevel = 0.5;

Seq.tPulseRotate = HW.tFlip180Def*2;
Seq.AmpPulseRotate = HW.TX.AmpDef;

% Seq.tAQfOffset = 10e-3;                % duration of frequency tracking (pay attention to chemical shift oscillation)
% Seq.nPeriods = 1;                      % number of periods
% Seq.SearchFlipDegStepsPerPeriod = 12;  % number of AQ windows per period

Seq.excitationPulse = @Pulse_Rect;  % pulse shape function of 90 degree excitation pulse

[HW, mySave, Seq] = Find_ExcitationPulse(HW, mySave, Seq);


%% CPMG echo train for 180 degrees pulse amplitude and phase
LoadSystem;                           % load system parameters
HW.FindFrequencyPause = T1*2;
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);

HW.FindFrequencyPause = max(500e-3, T1*2);

Seq.T1 = max(500e-3, T1*2);
Seq.tExcitation = 100e-6;
Seq.solidEchoPulse = @Pulse_Rect;  % other pulse function for no solid echo
% Seq.solidEchoPulse = @Pulse_Rect_SolidEcho;  % pulse shape function of solid echo pulse
% Seq.tauSolidEcho = 10e-6;
Seq.tEcho = 1e-3;
% Seq.refocusingPulse = @Pulse_Rect;
Seq.tRefocus = 200e-6;  %*Seq.refocusingPulse(HW, 'Time');
Seq.tAQEcho = min(500e-6, Seq.tEcho-Seq.tRefocus-100e-6);  % Seq.tRefocus * pi/2;
Seq.tEchoTrain = min(100e-3, T1*2);

Seq.Find_Frequency_interval = 0;
% find frequency more exactly with smaller bandwidth
HW.FindFrequencySweep.fOffsetFIDsStdMaxValue = 20;
HW.FindFrequencySweep.span = 0e3;
HW.FindFrequencySweep.nMeasurements = 1;
HW.FindFrequencySweep.tPulseDivider = 1;
HW.FindFrequencySweep.fSample = 2e3;
HW.FindFrequencySweep.nSamples = 40;
HW.FindFrequencySweep.maxTime = 0;
HW.FindFrequencySweep.doPlot = 1;


[HW, mySave, Seq] = Find_RefocusingPulse(HW, mySave, Seq);

LoadSystem;                           % load new parameters from file.

% search frequency with new values
[HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);


%% -----------------------------------------------------------------------------
% (C) Copyright 2024-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
