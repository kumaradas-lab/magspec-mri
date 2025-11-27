function [data, SeqOut, mySave] = sequence_RecoveryCPMG(HW, Seq, mySave)
%% Inversion or saturation followed by a CPMG train to measure T1 vs. T2
%
%   [data, SeqOut, mySave] = sequence_RecoveryCPMG(HW, Seq, mySave)
%
% This sequence acquires CPMG Echo trains optionally preceeded by preparation
% pulses.  It can be used to measure T1, T2 or T1-T2 correlation maps in 0D or
% in a slice.
%
% INPUT:
%   HW      HW structure or object.
%
%   Seq     A structure with the following optional fields. If the fields are
%           omitted or empty, default values are used:
%     inversionPulse
%             Handle to a pulse shape function for the 180 degrees inversion
%             pulse before the CPMG train (default: @Pulse_Rect).
%     saturationPulse
%             Handle to a pulse shape function for the 90 degrees saturation
%             pulse before the CPMG train (default: @Pulse_Rect).
%     excitationPulse
%             Handle to a pulse shape function for the 90 degrees excitation
%             pulse for the CPMG (default: @Pulse_Rect).
%     refocusingPulse
%             Handle to a pulse shape function for the 180 degrees refocusing
%             pulse for the CPMG (default: @Pulse_Rect).
%     excitationAngle
%             excitation angle of the (90 degrees) excitation pulse in degrees
%             (default: 90)
%     refocusingAngle
%             excitation angle of the (180 degrees) refocusing pulses in degrees
%             (Default: HW.RecoveryCPMG.refocusingAngle or 180)
%     refocusingPhase
%             Phase of the refocusing pulse (relative to the excitation pulse)
%             in degrees.
%             (Default: HW.RecoveryCPMG.refocusingPhase or 90)
%     refocusingPhaseOffset
%             Additional phase offset of the refocusing pulses in degrees. This
%             can be:
%               * A scalar that applies to all refocusing pulses in the echo
%                 train.
%               * A vector of size 1xSeq.nEcho.
%               * A string "alternateIncrement" in which case the refosing
%                 direction alternates after each second pulse.
%             (Default: 0)
%     tExcitation
%             duration of the excitation pulse in s. If this is set, the
%             amplitude of the excitation pulse is set such that the pulse flips
%             "excitationAngle" degrees. Otherwise, the amplitude is chosen such
%             that the refocusing pulse of the same duration flips
%             "refocusingAngle" degrees. (default: [])
%     tRefocus
%             duration of the refocusing pulse in s. If this is set, the
%             amplitude of the refocusing pulse is set such that the pulse flips
%             "refocusingAngle" degrees. Otherwise, the default amplitude is
%             used and the duration is set accordingly.
%             (default: Seq.tExcitation)
%     tPreheatRfAmp
%             Turn on the rf amplifier just before the excitation (or
%             preparation) pulse for this time in seconds.
%             (default: 0)
%     AQPhaseOffset
%             phase offset of the receiver for each Echo train in degrees
%             (default: 0)
%     TXPhaseOffset
%             phase offset of the refocusing pulses of each Echo train in
%             degrees (default: 0)
%
%     Find_Frequency_interval
%             Minimum time between two frequency sweeps in seconds. For this to
%             work properly, the function must be used with the optional input
%             and output argument mySave. Set to "Inf" if the frequency sweep
%             should be omitted. For further settings, see HW.FindFrequencySweep
%             (default: HW.FindFrequencySweep.maxTime if called with input
%             argument "mySave", Inf otherwise).
%
%     T2 (CPMG) SETTINGS:
%     T2Estimated
%             estimated T2 time of the sample in s
%             (default: 150e-3)
%     T2EstimatedMin
%             estimated minimum value of the range of T2 times of the sample
%             in s (default: Seq.T2Estimated/2)
%     T2EstimatedMax
%             estimated maximum value of the range of T2 times of the sample
%             in s (default: Seq.T2Estimated*2)
%     tEcho
%             Echo time of the CPMG Echo train in s. The set value is increased
%             to HW.RecoveryCPMG.tEchoMin if necessary.
%             (default: Seq.T2EstimatedMin/2 but at least 8 times the refocusing
%             pulse length)
%     tEchoTrain
%             (approximate) total time of Echo train in s. nEcho takes
%             precedence.
%             (default: Seq.T2EstimatedMax*2)
%     nEcho
%             number of Echoes in CPMG Echo train in s
%             (default: round(Seq.tEchoTrain/Seq.tEcho))
%     FixfLarmorTotEcho
%             Boolean value. If it evaluates to true, the frequency of the TX
%             pulses and the acquisition windows is adapted such that the Echo
%             time is an integer multiple of its period. This should only be
%             used in conjunction with slice selection. (Default: false)
%             The actually used frequency is returned in SeqOut.fLarmor.
%             See also "get_fLarmorFitTotEcho".
%     SteadyState_PreShots180
%             Number of Echoes without acquisition window at beginning of each
%             CPMG train.
%     tAQEcho
%             acquisition window duration in seconds of Echoes in CPMG train.
%             Set to 0 for a single sample at the Echo times.
%             (default: max(1e-6, min(100e-6,
%             Seq.tEcho-2*max(get_DeadTimeTX2RX(HW,Seq.fSample),get_DeadTimeRX2TX(HW,Seq.fSample))-Seq.tAQEchoDelay-Seq.tRefocus))
%             or the value in HW.RecoveryCPMG.tAQEcho)
%     tAQEchoDelay
%             delays the acquisition windows in s
%             (default: HW.RX.RiseTime or the value in
%             HW.RecoveryCPMG.tAQEchoDelay)
%     fSample
%             Sample rate of acquisition windows at Hahn echoes in Hertz
%             (default: 500e3 or the value in HW.RecoveryCPMG.fSample)
%     tAQFID
%             Duration of the acquisition window after the excitation pulse in
%             seconds. This value must be smaller or equal to Seq.tEcho/4. If
%             the time is set to 0, a single sample of the FID is acquired. If
%             the time is set to Inf, the maximum possible time (until the first
%             inversion pulse of the CPMG echo train) is acquired. For negative
%             values, the acquisition of the FID is omitted. (default: -1)
%     tAQFIDStart
%             Time in seconds after the center of the excitation pulse when the
%             acquisition window of the FID should start. If this time is too
%             short, it is automatically extended to the lowest possible value.
%             (default: 0)
%     fSampleFID
%             Sample rate of acquisition windows at FID (or solid echo) in Hertz
%             (default: Seq.fSample or the value in HW.RecoveryCPMG.fSampleFID)
%     AQFIDMean
%             Boolean value. If true, the average of the acquired FID is used
%             for the T2 fit. Otherwise, each sample of the FID is used
%             (default: false).
%     tauSolidEcho
%             Time between first (excitation) and second ("refocusing") 90
%             degree pulses for solid echo in seconds. See Pulse_Rect_SolidEcho.
%             (Default: 2*HW.tFlip90Def)
%
%     T1 SETTINGS:
%     Recovery
%             The type of the recovery experiment. Either 'Inversion' (180
%             degree preparation pulse), 'Saturation' (90 degree preparation
%             pulse) or 'Decay' (e.g. Spin-Lock preparation).
%             (default: 'Inversion')
%     tRelax
%             relaxation time after the CPMG Echo train in s
%             (default: Seq.T1EstimatedMax*5)
%     tSlice
%             time between CPMG echo trains of different slices in s
%             (default: Seq.tRelax)
%     T1Estimated
%             estimated T1 time of the sample in s
%             (default: 200e-3)
%     T1EstimatedMin
%             estimated minimum value of the range of T1 times of the sample
%             in s (default: Seq.T1Estimated/3)
%     T1EstimatedMax
%             estimated maximum value of the range of T1 times of the sample
%             in s (default: Seq.T1Estimated*3)
%     nTau1SteadyState
%             number of preparation-CPMG train sequences without data
%             acquisition before the actual measurement (default: 1)
%     Tau1Start
%             shortest recovery time for automatic spacing in s
%             (default: Seq.T1EstimatedMin/3)
%     Tau1End
%             longest recovery time for automatic spacing in s
%             (default: Seq.T1EstimatedMax*3)
%     Tau1Log
%             Boolean: If true, the nTau1 recovery times between Tau1Start and
%             Tau1End are spaced logarithmically. They are linearly spaced
%             otherwise.
%     Tau1StartFactor
%             increment factor for steps in recovery times. (default: 1 if
%             Seq.Tau1Log is false or exp(1) if Seq.Tau1Log is true)
%     nTau1
%             number of different recovery times. Set to 0 to acquire a CPMG
%             Echo train without prior preparation.
%             (default: round(Seq.Tau1End / (Seq.Tau1Start*Seq.Tau1StartFactor))
%             if Seq.Tau1Log is false or
%             round(1 + (log10(Seq.Tau1End) - log10(Seq.Tau1Start)) / ...
%                       (log10(Seq.Tau1Start*Seq.Tau1StartFactor)-log10(Seq.Tau1Start))))
%             if Seq.Tau1Log is true. But at least 3 in either case.
%     Tau1
%             Array of recovery times (time between center of preparation pulse
%             and center of excitation pulse) in s. If this is omitted or empty
%             and nTau1 > 0, the following defaults are used:
%             linspace(Seq.Tau1Start, Seq.Tau1End, Seq.nTau1)
%             if Seq.Tau1Log is false or
%             logspace(log10(Seq.Tau1Start), log10(Seq.Tau1End), Seq.nTau1)
%             if Seq.Tau1Log is true.
%     tShiftPreparation
%             Additional time between center of preparation pulse and center of
%             excitation pulse in seconds. This time is *not* included in the
%             evaluation of T1 which uses Tau1. This time can be used to account
%             for durations during the preparation time during which the
%             relevant magnetization doesn't change (e.g., while the spin-system
%             is *not* aligned with B0 during the preparation time).
%             (Default: 0)
%     PulsePreparation
%             Structure with settings for the pulse shape function of the
%             preparation pulses. (See Seq.saturationPulse or Seq.inversionPulse
%             above)
%     tSaturation
%             Duration of the saturation pulse in seconds. By default, the
%             duration is such that its bandwidth matches the excitation pulse.
%             This value is only used if Seq.Recovery is 'Saturation' or
%             'Decay'.
%     SpoilSaturation
%             In case Seq.Recovery is set to 'Saturation', a spoiler gradient
%             pulse can be used to reduce the influence of the excitation effect
%             of the preparation pulse. This Boolean value controls whether a
%             spoiler gradient pulse is used after the preparation pulse.
%             (Default: true)
%     tSaturationGrad
%             Duration of the (optional) spoiler gradient pulse after the
%             saturation pulse in s (default: 3e-3)
%     AmpSaturationGrad
%             Amplitude of the (optional) spoiler gradient pulse after the
%             saturation pulse in T
%             (default: min(HW.Grad.MaxAmp(1:3))/2)
%     GetDataEarly
%             If true, the data is supplied by the MRT early (packet end), i.e.
%             at the end of the last acquisition window (of each average).
%             (Default: Seq.tRelax > 1)
%     Function_Prepare_Measurement
%             A function handle to a prepare function with the following
%             signature:
%                   [HW, Seq, AQ, TX, Grad] = @(HW, Seq, AQ, TX, Grad)
%             The function is called after the basic pulse program (without
%             averages or phase cycling) is created. All elements have the
%             sorting nActions x tReps x nTau1.
%             It can be used to add more sophisticated preparations.
%             Set Seq.T2Prepare.tRepsPrepare accordingly if the number of tReps
%             used in the preparation is other than 1. Adjust Seq.tRep,
%             Seq.tOffset, Seq.CLTime, Seq.DigitalIO, TX, AQ, and Grad
%             accordingly.
%
%     EVALUATION SETTINGS:
%     Plot
%             Boolean: plot the acquired data on a linear (non-wrapped) time
%             axis (default: false)
%     PlotTR
%             Boolean: plot the acquired data where all tReps are stacked upon
%             each other (default: false)
%     FitT1
%             Boolean: fit T1 time to values of 1st Echo in each recovery
%             (default: Seq.nTau1>2)
%     PlotT1
%             Figure number for the T1 plot. If set to 0, no figure is shown.
%             (Default: 101)
%     FitT2
%             Boolean: fit T2 time of CPMG with FitT2AtTau1 (see below)
%             (default: Seq.nEcho>3)
%     FitT2AtTau1
%             indices for the recovery times that are selected for T2 time fit
%             (default: max(1,Seq.nTau1))
%     PlotT2
%             Figure number for the (first) T2 plot. If set to 0, no figure is
%             shown. (Default: 90)
%     PlotT1T2
%             Figure number for the (first) plot with T1-T2-map data. If set to
%             0, no figure is shown. (Default: 80)
%     SubtractEmptyMeasurement
%             Boolean: If true, the Echo amplitudes of an empty measurement are
%             subtracted from the actual measurement. (default: false)
%     PostProcessSequenceLocal
%             If false, the following and the fullScaleReference corrections are
%             omitted. (default: Seq.PostProcessSequence)
%     fitExpCorr
%             Structure with settings for fit_exp that is used for correcting
%             the data. See function "fit_exp" for available settings.
%             (Default: correct only frequency and phase offset if not using
%             slice, no fitting, no plots)
%     fitExpT1
%             Structure with settings for fit_exp that is used for fitting the
%             T1 data. See function "fit_exp" for available settings. (Default:
%             don't re-run the corrections, use "EndOffset", use settings from
%             fitExpCorr for the rest)
%     fitExpT2
%             Structure with settings for fit_exp that is used for fitting the
%             CPMG Echo train(s). See function "fit_exp" for available settings.
%             (Default: don't re-run the corrections, don't use "EndOffset", use
%             settings from fitExpCorr for the rest)
%
%     SLICE GRADIENT SETTINGS:
%     thicknessSlice
%             Thickness of slice selected by continuous gradient in m. The
%             slice gradient amplitude is calculated with respect to the 180
%             degree refocusing pulse length (the effective slice profile of the
%             shorter 90 degree excitation pulse is wider). If set to Inf, no
%             slice gradient is used. (default: Inf)
%     MaxGradAmpSlice
%             Maximum amplitude of the slice gradient in T/m. This leads to
%             potentially longer rf pulses. (default: HW.Grad.MaxAmpSlice)
%     Slice.offBetweenCPMG
%             Boolean: If the relaxation pause between the CPMG Echo trains long
%             enough, the slice gradient is turned off. (default: true)
%     Slice.offAfterPreparation
%             Boolean: For recovery experiments, if the time between preparation
%             pulse (inversion or excitation pulse) and the CPMG Echo train is
%             long enough, the slice gradient is turned off. (default: false)
%     Slice.nRamp
%             Number of steps in ramp. (default: 1)
%
%     AVERAGING SETTINGS:
%     SeqAverage
%             Structure with the following (optional) fields:
%       average
%               Number of times that the measurement for each preparation time
%               is repeated. (default: 1)
%       averageBreak
%               Additional time in seconds after the block with a given
%               preparation time was measured. (default: 0)
%       SaveSeqAverageData
%               Boolean value to indicate whether to sace data of all averages.
%               (default: true)
%       TXPhasePrepareIncrement
%               Phase increment in degrees of preparation pulses (inversion or
%               saturation). (default: [0, 0, 180, 0] if Seq.SeqAverage.average
%               is a multiple of 4, 0 otherwise)
%       TXPhaseExcitationIncrement
%               Phase increment of excitation pulses (of CPMG) from one average
%               to the next. (default: 180 for phase cycling)
%       TXPhaseRefocusIncrement
%               Phase increment of refocusing pulses (of CPMG) from one average
%               to the next. (default: 0)
%       AQPhaseIncrement
%               Phase increment of acquisitions (at echoes of CPMG) from one
%               average to the next.
%               (default: Seq.SeqAverage.TXPhaseExcitationIncrement)
%       RandomTXRXPhaseOffset
%               Boolean value to indicate whether a random phase should be added
%               to all rf pulses and acquisitions from one average to the next.
%               (default: false)
%       TXPhasePrepareOffset
%               Additional phase offset of preparation pulses in degrees.
%               (default: 0)
%       TXPhaseExcitationOffset
%               Additional phase offset of excitation pulses (of CPMG) in
%               degrees. (default: 0)
%       TXPhaseRefocusOffset
%               Additional phase offset of refocusing pulses (of CPMG) in
%               degrees. (default: 0)
%       AQPhaseOffset
%               Additional phase offset of acquisitions (at echoes of CPMG) in
%               degrees. (default: 0)
%       GetDataAfterAverages
%               Collect the data that has been acquired so far immediately after
%               the averages for a given preparation time are done.
%
%
%     SAMPLE LIFT SETTINGS:
%     Lift
%             Structure with settings for a sample lift.
%       useLift
%               Boolean: If true, add a trigger signal to the end of the
%               sequence to start a pre-defined lift movement. If false, the
%               values of the other fields are ignored. (default: false)
%       skipLift
%               Boolean: Tread sequence like the lift is used but don't actually
%               send the Digital IO signal. (default: false)
%       startSetChannel
%               Channel number of the IO port that is connected to the "start
%               set" channel of the sample lift controller. (default: 3)
%       setBit0Channel
%               Channel number of the IO port that is connected to the channel
%               that sets the bit 0 of the set selector on the sample lift
%               controller. (default: 4)
%       useSet
%               Number of the saved set that is triggered - either 1 or 2.
%               (default: 1)
%       checkEncoder
%               Check whether motor position and encoder position coincide.
%               This errors if the sample lift is still moving when its position
%               is querried. (default: true)
%
%
%     COIL DAMPING SETTINGS:
%     DampCoil
%             Structure for setting the coil damping circuit. These settings are
%             deprecated. Consider using HW.TX.DampCoil instead. The following
%             fields are used:
%       Enable
%               Boolean: If true, damp coil after 180 degree refocusing pulses.
%               (default: false)
%       Delay
%               Delay after 180 degree refocusing pulses before damping signal
%               in s (default: 0)
%       Duration
%               Duration of damping signal after 180 degree refocusing pulses
%               in s (default: 9e-6)
%       DOChannel
%               Number of the Digital Output channel that switches damping
%               circuit. (default: 1)
%       DampExcitation
%               Boolean: If true, damp coil after 90 degree excitation pulses.
%               (default: false)
%       DelayExcitation
%               Delay after 90 degree excitation pulses before damping signal
%               in s (default: 0)
%       DurationExcitation
%               Duration of damping signal after 90 degree excitation pulses
%               in s (default: 9e-6)
%
%
%     EXTERNAL POWER SUPPLY SETTINGS:
%     PowerSupply
%             Structure with settings for an external power supply that is used
%             to power the slice gradient. The following fields can be used:
%       EnableAnalogChannel
%               Number of the digital output channel that is used to set the
%               external power supply to analog remote control. If set to 0, no
%               signal is emitted.
%               (Default: HW.PowerSupply.EnableAnalogChannel if available, 0
%               otherwise)
%       EnableAnalogTOffset
%               Offset time in seconds before ramping up the slice gradient when
%               the enable analog control signal is switched on.
%               (Default: HW.PowerSupply.EnableAnalogTOffset if available, 0
%               otherwise)
%       EnableAnalogTPostset
%               Postset time in seconds after ramping down the slice gradient
%               when the enable analog control signal is switched off.
%               (Default: HW.PowerSupply.EnableAnalogTPostset if available, 0
%               otherwise)
%
%
%     FULL SCALE REFERENCE SETTINGS:
%     MeasureFullScaleReference
%             Use measurement as a reference for full scale
%             (default: 0)
%     fullScaleReference.Amplitude
%               amplitude that was measured for a 100% water sample
%               (default: [])
%     fullScaleReference.sampleCrossSection
%               cross section of the reference water sample in m^2
%               (default: [])
%     fullScaleReference.thicknessSlice
%               slice thickness that was used in the reference measurement in m
%               (default: [])
%     sampleCrossSection
%             cross section of the actual sample in m^2
%             (default: [])
%
%   mySave
%           structure as created by LoadSystem. This argument is optional. But
%           it is necessary to be passed when working with
%           Seq.Find_Frequency_interval.
%
% OUTPUT:
%   data
%           structure containing the measurement data and evaluation results.
%
%   SeqOut
%           Same as input structure "Seq" but with the following added fields:
%
%       amplitude2fullScaleReference
%             Factor that was used to convert amplitude to water parts.
%
%   mySave
%           structure with updated fields
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012,2015-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


% FIXME: Use add_DigitalIO


%% default parameters
if nargin < 3, mySave = struct(); end

if isemptyfield(Seq, 'PreProcessSequence'), Seq.PreProcessSequence  = 1;  end
if isemptyfield(Seq, 'StartSequence'),      Seq.StartSequence       = 1;  end
if isemptyfield(Seq, 'PollPPGfast'),        Seq.PollPPGfast         = 1;  end
if isemptyfield(Seq, 'GetRawData'),         Seq.GetRawData          = 1;  end
if isemptyfield(Seq, 'PostProcessSequence'),Seq.PostProcessSequence = 1;  end
if isemptyfield(Seq, 'PostProcessSequenceLocal'), Seq.PostProcessSequenceLocal = Seq.PostProcessSequence; end


if isemptyfield(Seq, 'Find_Frequency_interval')
  if nargin < 3
    Seq.Find_Frequency_interval = Inf;
  else
    Seq.Find_Frequency_interval = HW.FindFrequencySweep.maxTime;
  end
end


if Seq.StartSequence
  %% Frequency Sweep
  if isemptyfield(mySave, 'lastTime'), mySave.lastTime = 0; end
  if isinf(Seq.Find_Frequency_interval) || ...
      (now*24*3600-mySave.lastTime < Seq.Find_Frequency_interval-0.01)
    % Use last results of frequency sweep (i.e. mySave.HW.B0)
    [HW, mySave] = Find_Frequency_Sweep(HW, mySave, Inf);
  else
    [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);  % Find magnet frequency
  end
end


if Seq.PreProcessSequence
  if isemptyfield(Seq, 'PlotTR'),         Seq.PlotTR        = 0;                  end
  if isemptyfield(Seq, 'Plot'),           Seq.Plot          = 0;                  end

  if isemptyfield(Seq, 'thicknessSlice'),  Seq.thicknessSlice = Inf;  end  % thickness of slice selected by continuous gradient in m
  iDevice = 1;  % FIXME: Support pulse program at secondary devices.
  if isemptyfield(Seq, 'MaxGradAmpSlice')
    % maximum amplitude for slice gradient
    Seq.MaxGradAmpSlice = HW.Grad(iDevice).MaxAmpSlice;
  end

  if isemptyfield(Seq, 'excitationPulse'), Seq.excitationPulse = @Pulse_Rect_Slice; end % pulse shape function of 90 degree excitation pulse of CPMG
  if isemptyfield(Seq, 'inversionPulse'),  Seq.inversionPulse  = @Pulse_Rect_Slice; end % pulse shape function of inversion pulse before CPMG
  if isemptyfield(Seq, 'saturationPulse'), Seq.saturationPulse = @Pulse_Rect_Slice; end % pulse shape function of saturation pulse before
  if isemptyfield(Seq, 'excitationAngle'), Seq.excitationAngle = 90;              end % excitation angle of the excitation pulses in degrees
  if ~isfield(Seq, 'tExcitation')
    % duration of the excitation pulses in seconds
    Seq.tExcitation = [];
  end
  if isemptyfield(Seq, 'tPreheatRfAmp'),  Seq.tPreheatRfAmp = 0;                  end % turn on the rf amplifier before the sequence for this time in s

  if isemptyfield(Seq, 'SpoilSaturation'),  Seq.SpoilSaturation = true;  end  % use spoiler after saturation pulse
  if isemptyfield(Seq, 'tSaturationGrad'),  Seq.tSaturationGrad = 3e-3;  end
  if isemptyfield(Seq, 'AmpSaturationGrad'),  Seq.AmpSaturationGrad = min(HW.Grad(iDevice).MaxAmp(1:3))/2;  end

  if isemptyfield(Seq, 'AQPhaseOffset'),  Seq.AQPhaseOffset = 0;                  end % phase offset of the receiver for each Echo train
  if isemptyfield(Seq, 'TXPhaseOffset'),  Seq.TXPhaseOffset = 0;                  end % phase offset of the refocusing pulses (180er of CPMG) of each Echo train

  % T2 (CPMG) settings
  if isemptyfield(Seq, 'T2Estimated'),    Seq.T2Estimated   = 150e-3;             end % estimated T2 time of sample in s
  if isemptyfield(Seq, 'T2EstimatedMin'), Seq.T2EstimatedMin = Seq.T2Estimated/2; end % estimated minimum value of T2 values of sample in s
  if isemptyfield(Seq, 'T2EstimatedMax'), Seq.T2EstimatedMax = Seq.T2Estimated*2; end % estimated maximum value of T2 values of sample in s
  if isemptyfield(Seq, 'refocusingPulse')
    % pulse shape function of 180 degree refocusing pulse of CPMG
    Seq.refocusingPulse = @Pulse_Rect_Slice;
  end
  if isemptyfield(Seq, 'refocusingAngle')
    % excitation angle of the refocusing pulses in degrees
    if isemptyfield(HW, {'RecoveryCPMG', 'refocusingAngle'})
      Seq.refocusingAngle = 180;
    else
      Seq.refocusingAngle = HW.RecoveryCPMG.refocusingAngle;
    end
  end
  if isemptyfield(Seq, 'tRefocus')
    % duration of the refocusing pulses in seconds
    Seq.tRefocus = Seq.tExcitation;
  end
  if isempty(Seq.tExcitation) && ~isemptyfield(HW, {'RecoveryCPMG', 'tFlip90Def'})
    Seq.tExcitation = HW.RecoveryCPMG.tFlip90Def;
  end
  if isempty(Seq.tRefocus)
    if ~isempty(HW.RecoveryCPMG.tFlip180Def)
      tFlip180Def = HW.RecoveryCPMG.tFlip180Def;
    else
      tFlip180Def = HW.tFlip180Def;
    end

    Seq.tRefocus = tFlip180Def * (Seq.refocusingAngle/180) * Seq.refocusingPulse(HW, 'Amp');

    if ~isinf(Seq.thicknessSlice)
      Seq.tRefocus = max(Seq.tRefocus, ...
        1/Seq.MaxGradAmpSlice * Seq.refocusingPulse(HW, 'Time') / (HW.GammaDef/(2*pi) * Seq.thicknessSlice));
    end
  end

  if isemptyfield(Seq, 'tEcho'),  Seq.tEcho = max(8*Seq.tRefocus, Seq.T2EstimatedMin/3); end  % Echo time in s
  if ~isemptyfield(HW, {'RecoveryCPMG', 'tEchoMin'}) && Seq.tEcho < HW.RecoveryCPMG.tEchoMin
    warning('PD:sequence_RecoveryCPMG:tEchoShortened', ...
      'Seq.tEcho (%.3f ms) has been increased to HW.RecoveryCPMG.tEchoMin (%.3f ms).', ...
      Seq.tEcho*1e3, HW.RecoveryCPMG.tEchoMin*1e3);
  end
  if isemptyfield(Seq, 'tEchoTrain'),     Seq.tEchoTrain    = Seq.T2EstimatedMax*3; end % total time of Echo train in s
  if isemptyfield(Seq, 'nEcho'),          Seq.nEcho         = round(Seq.tEchoTrain/Seq.tEcho);end % number of Echoes in CPMG

  if isemptyfield(Seq, 'refocusingPhase')
    if isemptyfield(HW, {'RecoveryCPMG', 'refocusingPhase'})
      Seq.refocusingPhase = 90;
    else
      Seq.refocusingPhase = HW.RecoveryCPMG.refocusingPhase;
    end
  end
  if isemptyfield(Seq, 'refocusingPhaseOffset')
    Seq.refocusingPhaseOffset = 0;
  end
  if ischar(Seq.refocusingPhaseOffset)
    switch Seq.refocusingPhaseOffset
      case 'alternateIncrement'
        Seq.refocusingPhaseOffset = mod(cumsum(90-90*(-1).^(1:Seq.nEcho)), 360);
      otherwise
        error('PD:sequence_RecoveryCPMG', ...
          'Setting Seq.refocusingPhaseOffset to "%s" is not supported.', ...
          Seq.refocusingPhaseOffset);
    end
  end
  if isemptyfield(Seq, 'FixfLarmorTotEcho'), Seq.FixfLarmorTotEcho = false;       end % adjust the Larmor frequency to the Echo time
  if Seq.thicknessSlice > 0.1 && Seq.FixfLarmorTotEcho
    warning('PD:sequence_RecoveryCPMG:NoSliceWithFixfLarmor', ...
      '"Seq.FixfLarmorTotEcho" should only be true if used in conjunction with a slice selection.');
  end
  if isemptyfield(Seq, 'SteadyState_PreShots180'), Seq.SteadyState_PreShots180 = 0; end % number of Echoes without AQ window at beginning of CPMG train
  if isemptyfield(Seq, 'fSample') % sample rate of acquisition windows at echoes
    if ~isemptyfield(HW.RecoveryCPMG, 'fSample')
      Seq.fSample = HW.RecoveryCPMG.fSample;
    else
      Seq.fSample = 500e3; % min(1000e3, max(20e3, 20/Seq.tAQEcho));
    end
  end
  Seq.fSample = 1/(round(HW.RX(iDevice).fSample/Seq.fSample)/HW.RX(iDevice).fSample);
  if isemptyfield(Seq, 'fSampleFID') % sample rate of acquisition windows at FID (or solid echo)
    if ~isemptyfield(HW.RecoveryCPMG, 'fSampleFID')
      Seq.fSampleFID = HW.RecoveryCPMG.fSampleFID;
    else
      Seq.fSampleFID = Seq.fSample;
    end
  end
  Seq.fSampleFID = 1/(round(HW.RX(iDevice).fSample/Seq.fSampleFID)/HW.RX(iDevice).fSample);

  if isemptyfield(Seq, 'tAQEchoDelay') % delays the acquisition windows in s
    if ~isemptyfield(HW.RecoveryCPMG, 'tAQEchoDelay')
      Seq.tAQEchoDelay = HW.RecoveryCPMG.tAQEchoDelay;
    elseif ~isempty(HW.RX(iDevice).RiseTime)
      Seq.tAQEchoDelay = HW.RX(iDevice).RiseTime;
    else
      Seq.tAQEchoDelay = 0;
    end
  end
  if isemptyfield(Seq, 'tAQEcho')  % acquisition window duration in seconds of Echoes in CPMG
    if ~isemptyfield(HW.RecoveryCPMG, 'tAQEcho')
      Seq.tAQEcho = HW.RecoveryCPMG.tAQEcho;
    else
      Seq.tAQEcho = max(1e-6, min(100e-6, ...
        Seq.tEcho - 2*max(get_DeadTimeTX2RX(HW, Seq.fSample), get_DeadTimeRX2TX(HW, Seq.fSample)) - Seq.tAQEchoDelay - Seq.tRefocus));
    end
  end
  if isemptyfield(Seq, 'tAQFID'),         Seq.tAQFID        = -1;                 end % acquisition time of FID in s
  if isemptyfield(Seq, 'tAQFIDStart'),    Seq.tAQFIDStart   = 0;                  end % start of FID acquisition in s
  if isemptyfield(Seq, 'AQFIDMean'),      Seq.AQFIDMean     = false;              end

  if isemptyfield(Seq, 'tauSolidEcho'),   Seq.tauSolidEcho  = 2*HW.tFlip90Def;    end

  if ~isfield(Seq, 'CLTime'),             Seq.CLTime        = [];                 end

  if isemptyfield(Seq, 'CalculateFFTOfData'), Seq.CalculateFFTOfData = 0; end
  if ~isinf(Seq.tAQFID) && Seq.tAQFID+Seq.tAQFIDStart > Seq.tEcho/2
    % FIXME: This limit is somewhat arbitrary.
    error('PD:sequence_RecoveryCPMG:tAQFID', ...
      ['Seq.tAQFID+Seq.tAQFIDStart must be <= Seq.tEcho/2 (%.3f ms). ', ...
      'Additionally, account for TX pulses and receiver dead times!'], ...
      Seq.tEcho/2*1e3);
  end

  % T1 settings
  if ~isemptyfield(Seq, 'Tau1')
    % Seq.Tau1 has priority over other settings
    Seq.Tau1Start = Seq.Tau1(1);
    Seq.Tau1End = Seq.Tau1(end);
    Seq.nTau1 = numel(Seq.Tau1);
  end
  if isemptyfield(Seq, 'Recovery'),       Seq.Recovery      = 'Inversion';        end % type of recovery experiment, 'Inversion' or 'Saturation'
  if isemptyfield(Seq, 'T1Estimated'),    Seq.T1Estimated   = 200e-3;             end % estimated T1 time of sample in s
  if isemptyfield(Seq, 'T1EstimatedMin'), Seq.T1EstimatedMin = Seq.T1Estimated/3; end % estimated minimum value of T1 values of sample in s
  if isemptyfield(Seq, 'T1EstimatedMax'), Seq.T1EstimatedMax = Seq.T1Estimated*3; end % estimated maximum value of T1 values of sample in s
  if isemptyfield(Seq, 'tRelax'),         Seq.tRelax        = Seq.T1EstimatedMax*5; end % relaxation time after CPMG Echo train in s
  if isemptyfield(Seq, 'tSlice'),         Seq.tSlice        = Seq.tRelax;         end % time between echo trains of different slices in s
  if isemptyfield(Seq, 'nTau1SteadyState'), Seq.nTau1SteadyState = 1;             end % number of pre-shots with 1st Tau1
  if isemptyfield(Seq, 'Tau1Start'),      Seq.Tau1Start     = Seq.T1EstimatedMin/3; end % shortest recovery time for automatic spacing in s
  if isemptyfield(Seq, 'Tau1End'),        Seq.Tau1End       = Seq.T1EstimatedMax*3; end % longest recovery time for automatic spacing in s
  if isemptyfield(Seq, 'Tau1Log'),        Seq.Tau1Log       = 1;                  end % Boolean to select logarithmic spacing of recovery times (linear spacing otherwise)
  if Seq.Tau1Log
    if isemptyfield(Seq, 'Tau1StartFactor'),  Seq.Tau1StartFactor = exp(1);       end % increment factor for steps
    if isemptyfield(Seq, 'nTau1') % number of recovery times. Set to 0 for CPMG only
      Seq.nTau1 = max(3, ...
        round(1 + (log10(Seq.Tau1End) - log10(Seq.Tau1Start)) / ...
                  (log10(Seq.Tau1Start*Seq.Tau1StartFactor)-log10(Seq.Tau1Start))));
    end
  else
    if isemptyfield(Seq, 'Tau1StartFactor'),  Seq.Tau1StartFactor = 1;            end % increment factor for steps
    if isemptyfield(Seq, 'nTau1') % number of recovery times. Set to 0 for CPMG only
      Seq.nTau1 = max(3, round(Seq.Tau1End / (Seq.Tau1Start*Seq.Tau1StartFactor)));
    end
  end
  if ~isfield(Seq, 'Tau1'),               Seq.Tau1          = [];                 end
  if isempty(Seq.Tau1) && (Seq.nTau1 > 0) % array of recovery times in s
    % create vector of recovery times
    if Seq.Tau1Log
      Seq.Tau1 = logspace(log10(Seq.Tau1Start), log10(Seq.Tau1End), Seq.nTau1);
    else
      Seq.Tau1 = linspace(Seq.Tau1Start, Seq.Tau1End, Seq.nTau1);
    end
  end

  if isemptyfield(Seq, 'tShiftPreparation')
    % Additional preparation time that doesn't "count" for the magnetization.
    % The same field name is used in `prepare_Recovery`.
    Seq.tShiftPreparation = 0;
  end

  if ~isfield(Seq, 'Function_Prepare_Measurement'),  Seq.Function_Prepare_Measurement = [];  end

  if isemptyfield(Seq, 'GetDataEarly'),  Seq.GetDataEarly   = Seq.tRelax > 1;     end % early data

  % sequence plot settings
  if ~isfield(Seq, 'plotSequence'),       Seq.plotSequence  = [];                 end
  % if isemptyfield(Seq.plotSequence, 'wraps'), Seq.plotSequence.wraps = Seq.nTau1SteadyState + Seq.SeqAverage.average; end
  if isemptyfield(Seq.plotSequence, 'wraps'), Seq.plotSequence.wraps = 1;         end
  if isemptyfield(Seq.plotSequence, 'xLim')
    if Seq.plotSequence.wraps > 1
      Seq.plotSequence.xLim = [-Inf, max([0 Seq.Tau1])+Seq.tEcho+Seq.tEchoTrain];
    end
  end

  % evaluation settings
  if isemptyfield(Seq, 'PlotT1'),         Seq.PlotT1        = 101;                end % figure number for the T1 plot
  if isemptyfield(Seq, 'PlotT2'),         Seq.PlotT2        = 90;                 end % first figure number for the T2 plot
  if isemptyfield(Seq, 'PlotT1T2'),       Seq.PlotT1T2      = 80;                 end % first figure number for the T1-T2 plots
  if isemptyfield(Seq, 'FitT1'),          Seq.FitT1         = Seq.nTau1>2;        end % Boolean: fit T1 time to values of 1st Echo in each recovery
  if isemptyfield(Seq, 'FitT2'),          Seq.FitT2         = Seq.nEcho>3;        end % Boolean: fit T2 time of CPMG with FitT2AtTau1 (see below)
  if isemptyfield(Seq, 'FitT2AtTau1')
    % indices for the recovery times that are selected for T2 time fit
    if strcmp(Seq.Recovery, 'Decay')
      Seq.FitT2AtTau1 = 1;
    else
      Seq.FitT2AtTau1 = max(1, Seq.nTau1);
    end
  end
  if isemptyfield(Seq, 'intensityProfileEchoNum'), Seq.intensityProfileEchoNum = 1; end % indices of the Echoes that are selected for the lift intensity profile
  if isemptyfield(Seq, 'intensityProfileTau1'), Seq.intensityProfileTau1 = max(1,Seq.nTau1); end % indices for the recovery times that are selected for the lift intensity profile
  if isemptyfield(Seq, 'SubtractEmptyMeasurement'), Seq.SubtractEmptyMeasurement = false; end % Boolean: If true, the Echo amplitudes of an empty measurement are subtracted from the actual measurement
  if isemptyfield(Seq, 'sampleCrossSection'), Seq.sampleCrossSection = 1;         end % sample Cross Section in m^2
  if isemptyfield(Seq, 'MeasureFullScaleReference'), Seq.MeasureFullScaleReference = 0; end
  if isemptyfield(Seq, 'fullScaleReference'), Seq.fullScaleReference = struct();  end

  if Seq.MeasureFullScaleReference
    Seq.fullScaleReference.sampleCrossSection = Seq.sampleCrossSection;
    Seq.fullScaleReference.thicknessSlice = Seq.thicknessSlice;
    Seq.fullScaleReference.tEcho = Seq.tEcho;
    Seq.fullScaleReference.nEcho = Seq.nEcho;
    Seq.fullScaleReference.tRelax = Seq.tRelax;
  elseif isfield(HW, 'RecoveryCPMG') && ~isemptyfield(HW.RecoveryCPMG, 'fullScaleReferencePath')
    [pathstr, name, ext] = fileparts(HW.RecoveryCPMG.fullScaleReferencePath);
    if exist(fullfile(pathstr, [name, HW.TX(iDevice).CoilName, ext]), 'file')
      fullScaleReference = load(fullfile(pathstr, [name, HW.TX(iDevice).CoilName, ext]), 'fullScaleReference');
      Seq.fullScaleReference = fullScaleReference.fullScaleReference;
      clear fullScaleReference
    end
  end

  if isemptyfield(Seq.fullScaleReference, 'SelectedTau1'),      Seq.fullScaleReference.SelectedTau1       = max(1,Seq.nTau1);           end
  if isemptyfield(Seq.fullScaleReference, 'thicknessSlice'),    Seq.fullScaleReference.thicknessSlice     = [];           end
  if isemptyfield(Seq.fullScaleReference, 'sampleCrossSection'),Seq.fullScaleReference.sampleCrossSection = [];           end
  if isemptyfield(Seq.fullScaleReference, 'Amplitude'),         Seq.fullScaleReference.Amplitude          = [];           end
  if isemptyfield(Seq.fullScaleReference, 'Factors'),           Seq.fullScaleReference.Factors            = [];           end
  if isemptyfield(Seq.fullScaleReference, 'nFactors'),          Seq.fullScaleReference.nFactors           = ceil(Seq.nEcho/10);           end
  if isemptyfield(Seq.fullScaleReference, 'PhaseOffset'),       Seq.fullScaleReference.PhaseOffset        = 0;           end
  if isemptyfield(Seq.fullScaleReference, 'fitStart'),          Seq.fullScaleReference.fitStart           = 5;           end
  if isemptyfield(Seq.fullScaleReference, 'hParent_fit_exp'),   Seq.fullScaleReference.hParent_fit_exp    = 1001;           end
  if isemptyfield(Seq.fullScaleReference, 'hParent_Factors'),   Seq.fullScaleReference.hParent_Factors    = 1002;           end

  if isempty(Seq.sampleCrossSection) || ...
      isempty(Seq.thicknessSlice) || ...
      isempty(Seq.fullScaleReference.thicknessSlice) || ...
      isempty(Seq.fullScaleReference.sampleCrossSection) || ...
      isempty(Seq.fullScaleReference.Amplitude) || ...
      Seq.MeasureFullScaleReference
    Seq.amplitude2fullScaleReference = 1;
  else
    Seq.amplitude2fullScaleReference = 1 / Seq.fullScaleReference.Amplitude ...
      * Seq.fullScaleReference.sampleCrossSection * min(Seq.fullScaleReference.thicknessSlice, abs(diff(HW.Grad(iDevice).ImageVol(3:4)))) ...
      / Seq.sampleCrossSection / min(Seq.thicknessSlice, abs(diff(HW.Grad(iDevice).ImageVol(3:4))));
  end

  if isemptyfield(Seq, 'fitExpCorr'),     Seq.fitExpCorr    = struct();           end
  fitExpCorr = Seq.fitExpCorr;
  if isemptyfield(fitExpCorr, 'hParent'),     fitExpCorr.hParent = 0;             end
  if isemptyfield(fitExpCorr, 'CorrectFrequencyDrift'), fitExpCorr.CorrectFrequencyDrift = 0; end
  if isemptyfield(fitExpCorr, 'CorrectFrequencyOffset')
    if isinf(Seq.thicknessSlice)
      fitExpCorr.CorrectFrequencyOffset = 1;
    else
      fitExpCorr.CorrectFrequencyOffset = 0;
    end
  end
  if isemptyfield(fitExpCorr, 'CorrectPhaseOffset')
    if isinf(Seq.thicknessSlice)
      fitExpCorr.CorrectPhaseOffset = 1;
    else
      fitExpCorr.CorrectPhaseOffset = 1;  % FIXME: Is this what we want?
    end
  end
  if isemptyfield(fitExpCorr, 'EndOffset')
    % FIXME: Improve this logic.
    if Seq.nTau1 > 1
      fitExpCorr.EndOffset = 1;
    else
      fitExpCorr.EndOffset = 0;
    end
  end
  if isemptyfield(fitExpCorr, 'RingFilter'),  fitExpCorr.RingFilter = 0;          end
  if isemptyfield(fitExpCorr, 'hasFID'),      fitExpCorr.hasFID = (Seq.tAQFID >= 0); end
  if isemptyfield(fitExpCorr, 'FIDMean'),     fitExpCorr.FIDMean = Seq.AQFIDMean; end
  if isemptyfield(fitExpCorr, 'SingleExp'),   fitExpCorr.SingleExp = 0;           end
  if isemptyfield(fitExpCorr, 'DoubleExp'),   fitExpCorr.DoubleExp = 0;           end
  Seq.fitExpCorr = fitExpCorr;

  if Seq.FitT1
    if isemptyfield(Seq, 'fitExpT1'),       Seq.fitExpT1      = struct();           end
    fitExpT1 = Seq.fitExpT1;
    if isemptyfield(fitExpT1, 'hParent'),     fitExpT1.hParent    = Seq.PlotT1;   end
    if isemptyfield(fitExpT1, 'CorrectFrequencyDrift'), fitExpT1.CorrectFrequencyDrift = 1; end
    if isemptyfield(fitExpT1, 'CorrectFrequencyOffset'), fitExpT1.CorrectFrequencyOffset = 0; end
    if isemptyfield(fitExpT1, 'CorrectPhaseOffset'), fitExpT1.CorrectPhaseOffset = 1; end
    if isemptyfield(fitExpT1, 'EndOffset')
      fitExpT1.EndOffset = ~strcmp(Seq.Recovery, 'Decay');
    end
    if isemptyfield(fitExpT1, 'hasFID'),      fitExpT1.hasFID     = 0;  end
    if isemptyfield(fitExpT1, 'RingFilter'),  fitExpT1.RingFilter = 0;  end
    if isemptyfield(fitExpT1, 'SingleExp'),   fitExpT1.SingleExp  = 1;  end
    if isemptyfield(fitExpT1, 'DoubleExp'),   fitExpT1.DoubleExp  = 1;  end
    Seq.fitExpT1 = fitExpT1;
  end

%   if Seq.FitT2
    if isemptyfield(Seq, 'fitExpT2'),       Seq.fitExpT2      = struct();           end
    fitExpT2 = Seq.fitExpT2;
    if isemptyfield(fitExpT2, 'hParent'),     fitExpT2.hParent  = Seq.PlotT1;     end
    if isemptyfield(fitExpT2, 'CorrectFrequencyDrift'), fitExpT2.CorrectFrequencyDrift = 0; end
    if isemptyfield(fitExpT2, 'CorrectFrequencyOffset'), fitExpT2.CorrectFrequencyOffset = 0; end
    if isemptyfield(fitExpT2, 'CorrectPhaseOffset'), fitExpT2.CorrectPhaseOffset = 0; end
    if isemptyfield(fitExpT2, 'EndOffset'),   fitExpT2.EndOffset = 0;             end
    if isemptyfield(fitExpT2, 'hasFID'),      fitExpT2.hasFID   = (Seq.tAQFID >= 0); end
    if isemptyfield(fitExpT2, 'FIDMean'),     fitExpT2.FIDMean  = 0;              end
    if isemptyfield(fitExpT2, 'RingFilter'),  fitExpT2.RingFilter = fitExpCorr.RingFilter; end
    if isemptyfield(fitExpT2, 'SingleExp'),   fitExpT2.SingleExp = 1;             end
    if isemptyfield(fitExpT2, 'DoubleExp'),   fitExpT2.DoubleExp = 1;             end
    Seq.fitExpT2 = fitExpT2;
%   end

  if isemptyfield(Seq, 'RawData'),        Seq.RawData       = ~or(Seq.Plot, Seq.PlotTR); end % Boolean: If true, run the evaluation on the un-corrected(?) raw data of the measurement

  % settings for slice gradient
  if ~isfield(Seq, 'Slice'), Seq.Slice = []; end
  if isemptyfield(Seq.Slice, 'offBetweenCPMG'), Seq.Slice.offBetweenCPMG = true;  end % Boolean: turn slice gradient off between CPMG Echo trains
  if isemptyfield(Seq.Slice, 'offAfterPreparation'), Seq.Slice.offAfterPreparation = false; end % Boolean: turn slice gradient off between preparation pulse and CPMG Echo train
  if isemptyfield(Seq.Slice, 'nRamp'),          Seq.Slice.nRamp = 1;              end % Number of steps in Ramp
  if Seq.Slice.nRamp < 1 || (Seq.Slice.nRamp ~= round(Seq.Slice.nRamp))
    error('PD:RecoveryCPMG:nRampNonInteger', '"Seq.Slice.nRamp" must be a positive integer.');
  end

  if HW.Grad(iDevice).Inductance(HW.Grad(iDevice).Slice.channel) == 0 && ...
      Seq.Slice.nRamp > 1
    if ~isinf(Seq.thicknessSlice)
      warning('PD:RecoveryCPMG:ZeroInductance', ...
        ['Number of segments for gradient ramps is >1 (%d). ', ...
        'But the inductance of the used gradient channel is 0. ', ...
        'Using a linear ramp with a single segment.'], Seq.Slice.nRamp);
    end
    Seq.Slice.nRamp = 1;
  end

  if ~isfield(Seq, 'Lift'), Seq.Lift = struct(); end
  if isemptyfield(Seq.Lift, 'useLift'),         Seq.Lift.useLift = false;         end % Boolean, move lift at end of sequence
  if isemptyfield(Seq.Lift, 'skipLift'),        Seq.Lift.skipLift = false;        end % Boolean, don't send digital IO signal for lift
  if isemptyfield(Seq.Lift, 'startSetChannel'), Seq.Lift.startSetChannel = 3;     end % digital output of drive-l that is connected to the start set port of the sample lift driver
  if isemptyfield(Seq.Lift, 'setBit0Channel'),  Seq.Lift.setBit0Channel = 4;      end % digital output of drive-l that is connected to the port that sets bit 0 of the set selection of the sample lift driver
  if isemptyfield(Seq.Lift, 'useSet'),          Seq.Lift.useSet = 1;              end % number of the saved set
  if isemptyfield(Seq.Lift, 'checkEncoder'),    Seq.Lift.checkEncoder = true;     end % check whether motor position and encoder position coincide

  if Seq.Lift.useLift && ...
      ((isa(HW, 'PD.HWClass') && ~isa(HW.Lift, 'PD.SampleLift')) || ...
      (~isa(HW, 'PD.HWClass') && (~isfield(HW, 'Lift') || ~isa(HW.Lift, 'PD.SampleLift'))))
    error('PD:sequence:NoSampleLift', 'Cannot use lift in sequence. It is not configured.');
  end

  % structure with settings for damping coil after TX pulses
  if isemptyfield(Seq, 'DampCoil'),  Seq.DampCoil = struct();  end
  if isemptyfield(Seq.DampCoil, 'Enable')
    % Damp coil after 180 degrees refocusing pulses
    Seq.DampCoil.Enable = false;
  end
  if isemptyfield(Seq.DampCoil, 'Delay')
    % Delay after 180 degrees refocusing pulse before damping signal in seconds
    Seq.DampCoil.Delay = 0;
  end
  if isemptyfield(Seq.DampCoil, 'Duration')
    % Duration of damp signal in seconds
    Seq.DampCoil.Duration = 9e-6;
  end
  if isemptyfield(Seq.DampCoil, 'DOChannel')
    % Digital Output channel that switches damping circuit
    Seq.DampCoil.DOChannel = 1;
  end
  if isemptyfield(Seq.DampCoil, 'DampExcitation')
    % Also damp after 90 degrees excitation pulse
    Seq.DampCoil.DampExcitation = Seq.DampCoil.Enable;
  end
  if isemptyfield(Seq.DampCoil, 'DelayExcitation')
    % Delay after 90 degrees excitation pulse before damping signal in seconds
    Seq.DampCoil.DelayExcitation = Seq.DampCoil.Delay;
  end
  if isemptyfield(Seq.DampCoil, 'DurationExcitation')
    % Duration of damp signal after 90 degrees excitation pulse
    Seq.DampCoil.DurationExcitation = Seq.DampCoil.Duration;
  end


  % structure with settings for external power supply
  if isemptyfield(Seq, 'PowerSupply'),  Seq.PowerSupply = struct(); end
  if isemptyfield(Seq.PowerSupply, 'EnableAnalogChannel')
    % digital output channel used for switching power supply to analog remote
    % control
    % 0 means signal is disabled
    if (isstruct(HW.PowerSupply) && isemptyfield(HW.PowerSupply, 'EnableAnalogChannel')) || ...
        isempty(HW.PowerSupply.EnableAnalogChannel)
      Seq.PowerSupply.EnableAnalogChannel = 0;
    else
      Seq.PowerSupply.EnableAnalogChannel = HW.PowerSupply.EnableAnalogChannel;
    end
  end
  if isemptyfield(Seq.PowerSupply, 'EnableAnalogTOffset')
    % offset of the analog remote control signal before the slice gradient pulse
    % in seconds
    if (isstruct(HW.PowerSupply) && isemptyfield(HW.PowerSupply, 'EnableAnalogTOffset')) || ...
        isempty(HW.PowerSupply.EnableAnalogTOffset)
      Seq.PowerSupply.EnableAnalogTOffset = 0;
    else
      Seq.PowerSupply.EnableAnalogTOffset = HW.PowerSupply.EnableAnalogTOffset;
    end
  end
  if isemptyfield(Seq.PowerSupply, 'EnableAnalogTPostset')
    % postset of the analog remote control signal after the slice gradient pulse
    % in seconds
    if (isstruct(HW.PowerSupply) && isemptyfield(HW.PowerSupply, 'EnableAnalogTPostset')) || ...
        isempty(HW.PowerSupply.EnableAnalogTPostset)
      Seq.PowerSupply.EnableAnalogTPostset = 0;
    else
      Seq.PowerSupply.EnableAnalogTPostset = HW.PowerSupply.EnableAnalogTPostset;
    end
  end


  if Seq.nTau1 > 0
    Seq.Tau1SteadyState = Seq.Tau1Start + zeros(1, Seq.nTau1SteadyState);
    Seq.Tau1All = [Seq.Tau1SteadyState, Seq.Tau1];
  else
    Seq.Tau1All = [];
  end


  Seq.tEchos = ((1:Seq.nEcho)*(Seq.tEcho*HW.MMRT(iDevice).fSystem))/HW.MMRT(iDevice).fSystem;  % time of echoes

  if ~isfield(Seq, 'SeqAverage'),         Seq.SeqAverage    = [];                 end % use average
  if isemptyfield(Seq.SeqAverage, 'average'),                   Seq.SeqAverage.average                  = 1;    end   % average number
  if isemptyfield(Seq.SeqAverage, 'averageBreak'),              Seq.SeqAverage.averageBreak             = 0;    end   % average break
  if isemptyfield(Seq.SeqAverage, 'SaveSeqAverageData'),        Seq.SeqAverage.SaveSeqAverageData       = 1;    end   % save data of all averages
  if isemptyfield(Seq.SeqAverage, 'TXPhasePrepareIncrement')
    % phase increment of preparation pulses (inversion or saturation)
    if mod(Seq.SeqAverage.average, 4) == 0
      Seq.SeqAverage.TXPhasePrepareIncrement = [0, 0, 180, 0];
    else
      Seq.SeqAverage.TXPhasePrepareIncrement = 0;
    end
  end
  if isemptyfield(Seq.SeqAverage, 'TXPhaseExcitationIncrement')
    % phase increment between averages of excitation pulse (90er of CPMG)
    Seq.SeqAverage.TXPhaseExcitationIncrement = 180;
  end
  if isemptyfield(Seq.SeqAverage, 'TXPhaseRefocusIncrement')
    % phase increment between averages of 180er of CPMG
    Seq.SeqAverage.TXPhaseRefocusIncrement = 0;
  end
  if isemptyfield(Seq.SeqAverage, 'AQPhaseIncrement')
    % phase increment between averages of acquisition window at Echoes
    Seq.SeqAverage.AQPhaseIncrement = Seq.SeqAverage.TXPhaseExcitationIncrement;
  end
  if isemptyfield(Seq.SeqAverage, 'RandomTXRXPhaseOffset')
    % Boolean: randomize TX pulse and acquisition phase (why?)
    Seq.SeqAverage.RandomTXRXPhaseOffset = 0;
  end
  if Seq.SeqAverage.RandomTXRXPhaseOffset
    [~, IX] = sort(rand(1, Seq.SeqAverage.average));
    Seq.SeqAverage.RandomTXRXPhaseOffset = linspace(0, 360-360./Seq.SeqAverage.average, Seq.SeqAverage.average);
    Seq.SeqAverage.RandomTXRXPhaseOffset = Seq.SeqAverage.RandomTXRXPhaseOffset(IX);
    Seq.SeqAverage.TXPhasePrepareOffset = Seq.SeqAverage.RandomTXRXPhaseOffset;
    Seq.SeqAverage.TXPhaseExcitationOffset = Seq.SeqAverage.RandomTXRXPhaseOffset;
    Seq.SeqAverage.TXPhaseRefocusOffset = Seq.SeqAverage.RandomTXRXPhaseOffset;
    Seq.SeqAverage.AQPhaseOffset = Seq.SeqAverage.RandomTXRXPhaseOffset;
  else
    if isemptyfield(Seq.SeqAverage, 'TXPhasePrepareOffset'),    Seq.SeqAverage.TXPhasePrepareOffset     = 0;  end % additional phase offset of preparation pulses
    if isemptyfield(Seq.SeqAverage, 'TXPhaseExcitationOffset'), Seq.SeqAverage.TXPhaseExcitationOffset  = 0;  end % additional phase offset of excitation pulse (90er of CPMG)
    if isemptyfield(Seq.SeqAverage, 'TXPhaseRefocusOffset'),    Seq.SeqAverage.TXPhaseRefocusOffset     = 0;  end % additional phase offset of 180er of CPMG
    if isemptyfield(Seq.SeqAverage, 'AQPhaseOffset'),           Seq.SeqAverage.AQPhaseOffset            = 0;  end % additional phase offset of acquisition windows
    % FIXME: Should the following conditions check the input size more specifically?
    if numel(Seq.SeqAverage.TXPhasePrepareOffset) < Seq.SeqAverage.average
      Seq.SeqAverage.TXPhasePrepareOffset = repmat(Seq.SeqAverage.TXPhasePrepareOffset, 1, Seq.SeqAverage.average/numel(Seq.SeqAverage.TXPhasePrepareOffset));
    end
    if numel(Seq.SeqAverage.TXPhaseRefocusOffset) < Seq.SeqAverage.average
      Seq.SeqAverage.TXPhaseRefocusOffset = repmat(Seq.SeqAverage.TXPhaseRefocusOffset, 1, Seq.SeqAverage.average/numel(Seq.SeqAverage.TXPhaseRefocusOffset));
    end
    if numel(Seq.SeqAverage.TXPhaseExcitationOffset) < Seq.SeqAverage.average
      Seq.SeqAverage.TXPhaseExcitationOffset = repmat(Seq.SeqAverage.TXPhaseExcitationOffset, 1, Seq.SeqAverage.average/numel(Seq.SeqAverage.TXPhaseExcitationOffset));
    end
    if numel(Seq.SeqAverage.AQPhaseOffset) < Seq.SeqAverage.average
      Seq.SeqAverage.AQPhaseOffset = repmat(Seq.SeqAverage.AQPhaseOffset, 1, Seq.SeqAverage.average/numel(Seq.SeqAverage.AQPhaseOffset));
    end
  end

  if isemptyfield(Seq.SeqAverage, 'GetDataAfterAverages')
    % collect data before starting the next preparation time
    Seq.SeqAverage.GetDataAfterAverages = false;
  end

  if isemptyfield(Seq, 'ResetDDSPhaseAtExcitation')
    % reset phase of transmission and acquisition at excitation pulse
    % (i.e. after preparation pulse)
    Seq.ResetDDSPhaseAtExcitation = 1;
  end


  % preliminary layout: [actions x CPMG x nTau1]
  Seq.tRep = repmat([zeros(1,Seq.nTau1>0), Seq.tEcho/2, repmat(Seq.tEcho,1,Seq.nEcho), Seq.tRelax], [1, 1, max(1,Seq.nTau1)+Seq.nTau1SteadyState]) + ...
    reshape([Seq.Tau1All(:).'+Seq.tShiftPreparation; zeros(1+Seq.nEcho+1, max(1,Seq.nTau1) + Seq.nTau1SteadyState)], [1, (Seq.nEcho+2+(Seq.nTau1>0)), (max(1,Seq.nTau1)+Seq.nTau1SteadyState)]);
  tRepsPerCMPG = size(Seq.tRep, 2);
  Seq.tOffset = zeros(size(Seq.tRep));

  [Grad(1:HW.Grad(iDevice).n).Time] = deal(NaN);
  [Grad(1:HW.Grad(iDevice).n).Amp] = deal(0);
  [Grad(1:HW.Grad(iDevice).n).Repeat] = deal([0, ones(1, length(Seq.tRep)-1)]);


  %% AQ
  if Seq.FixfLarmorTotEcho
    [Seq.fLarmor, Seq.tEcho, df] = get_fLarmorFitTotEcho(HW, HW.fLarmor, Seq.tEcho);
    if df > 10e3
      warning('PD:sequence_RecoveryCPMG:LargeDeviation:fLarmor', ...
        ['The used frequency (%.4f MHz) deviates significantly (%.2 kHz) ', ...
        'from the B0 Larmor frequency (%.4f MHz).'], ...
        Seq.fLarmor, df, HW.fLarmor);
    end
  else
    Seq.fLarmor = HW.fLarmor;
  end
  AQ.Frequency = Seq.fLarmor;


  %% TX (CPMG)
  % use excitation pulse with approximately the same bandwidth as the
  % refocusing pulses

  Seq.tRefocusBW = 1/(Seq.tRefocus) * Seq.refocusingPulse(HW, 'Time');
  if isempty(Seq.tExcitation)
    Seq.tExcitationBW = Seq.tRefocusBW / Seq.refocusingPulse(HW, 'Time') * Seq.excitationPulse(HW, 'Time');
    Seq.tExcitation = 1/Seq.tExcitationBW * Seq.excitationPulse(HW, 'Time');
  else
    Seq.tExcitationBW = 1/Seq.tExcitation * Seq.excitationPulse(HW, 'Time');
  end

  PulseExcitation.Bandwidth = Seq.tExcitationBW;
  PulseExcitation.FlipAngle = HW.tFlip90Def  * (Seq.excitationAngle/90*pi) / ...
    (HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef);
  PulseExcitation.MaxNumberOfSegments = 51;  % FIXME: Why?
  PulseExcitation.Frequency = AQ.Frequency(1);
  PulseExcitation.Phase = 0;
  % only applies to solid echo pulses
  PulseExcitation.tauSolid = Seq.tauSolidEcho;
  PulseExcitation.PhaseSolidInversionOffset = 90;
  PulseExcitation.iDevice = iDevice;
  Seq.pulseExcitation = Seq.excitationPulse(HW, 0, PulseExcitation);

  PulseRefocus.Bandwidth = Seq.tRefocusBW;
  PulseRefocus.FlipAngle = Seq.refocusingAngle/180*pi;
  PulseRefocus.MaxNumberOfSegments = 51;  % FIXME: Why?
  PulseRefocus.Frequency = AQ.Frequency(1);
  PulseRefocus.Phase = Seq.refocusingPhase;
  PulseRefocus.iDevice = iDevice;
  Seq.pulseRefocus = Seq.refocusingPulse(HW, 0, PulseRefocus);


  % preparation pulses
  if Seq.nTau1 > 0
    % preparation pulse
    if isemptyfield(Seq, 'PulsePreparation'), Seq.PulsePreparation = struct(); end
    if isemptyfield(Seq.PulsePreparation, 'MaxNumberOfSegments')
      Seq.PulsePreparation.MaxNumberOfSegments = 51;
    end
    if isemptyfield(Seq.PulsePreparation, 'MaxLength')
      Seq.PulsePreparation.MaxLength = Inf;
    end
    if isemptyfield(Seq.PulsePreparation, 'Frequency')
      Seq.PulsePreparation.Frequency = AQ.Frequency(1);
    end
    switch Seq.Recovery
      case 'Inversion'
        if isemptyfield(Seq, {'iLaplace2D', 'Recovery'})
          Seq.iLaplace2D.Recovery = 'Inversion';
        end
        Seq.tInvert = Seq.tRefocus / Seq.refocusingPulse(HW, 'Amp') * Seq.inversionPulse(HW, 'Amp');
        if isemptyfield(Seq.PulsePreparation, 'Bandwidth')
          Seq.PulsePreparation.Bandwidth = 1/Seq.tInvert * Seq.inversionPulse(HW, 'Time');
        end
        if isemptyfield(Seq.PulsePreparation, 'FlipAngle')
          % default to flip angle of refocussing pulses
          Seq.PulsePreparation.FlipAngle = Seq.refocusingAngle/180*pi;
        end
        if isemptyfield(Seq.PulsePreparation, 'Phase')
          Seq.PulsePreparation.Phase = 90;
        end
        Seq.pulsePrepare = Seq.inversionPulse(HW, 0, Seq.PulsePreparation);

      case {'Saturation', 'Decay'}
        if isemptyfield(Seq, {'iLaplace2D', 'Recovery'})
          Seq.iLaplace2D.Recovery = 'Saturation';
        end
        if isemptyfield(Seq, 'tSaturation')
          Seq.tSaturation = 1/Seq.tExcitationBW * Seq.saturationPulse(HW, 'Amp');
        end
        if isemptyfield(Seq.PulsePreparation, 'Bandwidth')
          Seq.PulsePreparation.Bandwidth = 1/Seq.tSaturation * Seq.saturationPulse(HW, 'Time');
        end
        if isemptyfield(Seq.PulsePreparation, 'FlipAngle')
          Seq.PulsePreparation.FlipAngle = pi/2;
        end
        if isemptyfield(Seq.PulsePreparation, 'Phase')
          Seq.PulsePreparation.Phase = 0;
        end
        % only applies to spin-locking pulses
        if isemptyfield(Seq.PulsePreparation, 'DurationSpinLock')
          Seq.PulsePreparation.DurationSpinLock = Seq.Tau1All;
        end
        Seq.pulsePrepare = Seq.saturationPulse(HW, 0, Seq.PulsePreparation);

        if Seq.SpoilSaturation
          % spoil saturation pulse
          if Seq.tSaturation/2 + HW.Grad(iDevice).tEC + Seq.tSaturationGrad + ...
              HW.Grad(iDevice).tRamp + HW.Grad(iDevice).tEC + Seq.tSaturation/2 + 50e-6 ...
              > Seq.Tau1Start
            error('increase Seq.Tau1Start or decrease Seq.tSaturationGrad');
          end
          for t=1:3
            if ~isinf(Seq.thicknessSlice) && t == HW.Grad(iDevice).Channel2xyzB(HW.Grad(iDevice).Slice.channel)
              continue;
            end
            Grad(t).Time = repmat(cumsum([Seq.tSaturation/2+HW.Grad(iDevice).tEC; HW.Grad(iDevice).tRamp; Seq.tSaturationGrad-HW.Grad(iDevice).tRamp*1; HW.Grad(iDevice).tRamp])*[1,nan(1,1+Seq.nEcho+1)], ...
              1, 1, Seq.nTau1 + Seq.nTau1SteadyState);
            Grad(t).Amp = repmat(Seq.AmpSaturationGrad*[0;1;1;0;]*[1,zeros(1,1+Seq.nEcho+1)], 1, 1, Seq.nTau1+Seq.nTau1SteadyState);
            Grad(t).Repeat = repmat([0, 0, ones(1, Seq.nEcho+1)], [1, 1, Seq.nTau1+Seq.nTau1SteadyState]);
          end
        end
      otherwise
        error('Please set Seq.Recovery to ''Inversion'', ''Saturation'', or ''Decay''!');
    end
  else
    Seq.pulsePrepare.Start = [];
  end


  if isempty(Seq.pulseExcitation.Start)
    % Use negative tOffset to allow longer preparation pulses
    Seq.tOffset(:,(Seq.nTau1>0)+1,:) = -60e-6;
  end


  % FIXME: Seq.TXPhaseOffset wird wohl in der eigentlichen Sequenz nicht
  % auf TX angewendet (only as a default for AQ).
  if isscalar(Seq.TXPhaseOffset)
    Seq.TXPhaseOffset = repmat(Seq.TXPhaseOffset,1,numel(Seq.tRep));
    TX.Repeat = repmat([zeros(1, Seq.nTau1>0), 0, zeros(1, Seq.nEcho>1), ones(1,max(0,Seq.nEcho-1)),0], [1, 1, max(1,Seq.nTau1)+Seq.nTau1SteadyState]);
  end
  if ischar(Seq.TXPhaseOffset)
    if strcmp(Seq.TXPhaseOffset, 'rand')
      Seq.TXPhaseOffset = repmat([zeros(1,Seq.nTau1>0),0,360*rand(1,Seq.nEcho),0], 1, 1, Seq.nTau1+Seq.nTau1SteadyState);
    end
    if strcmp(Seq.TXPhaseOffset, 'linear')
      Seq.TXPhaseOffset = repmat([zeros(1,Seq.nTau1>0),0,360*linspace(0,1,Seq.nEcho+1)], 1, 1, Seq.nTau1+Seq.nTau1SteadyState);
    end
    if strcmp(Seq.TXPhaseOffset, 'alternate')
      Seq.TXPhaseOffset = repmat([zeros(1,Seq.nTau1>0),0,0+90*(-1).^(1:Seq.nEcho),0], 1, 1, Seq.nTau1+Seq.nTau1SteadyState);
    end
  end
  % Seq.pulseRefocus.Phase = Seq.pulseRefocus.Phase + Seq.TXPhaseOffset;

  % initialize TX structure
  TX.Start = nan(max([size(Seq.pulsePrepare.Start,1)*(Seq.nTau1>0), length(Seq.pulseExcitation.Start) + length(Seq.pulseRefocus.Start)]), Seq.nEcho+2+(Seq.nTau1>0));
  TX.Duration = TX.Start;
  TX.Amplitude = TX.Start;
  TX.Frequency = TX.Start;
  TX.Phase = TX.Start;

  % Should be the frequency of the inversion pulses
  TX.Frequency(:) = AQ.Frequency(1);

  if Seq.nTau1 > 0
    % preparation pulse
    TX.Start(1:size(Seq.pulsePrepare.Start,1),1) = Seq.pulsePrepare.Start(:,1);
    TX.Duration(1:size(Seq.pulsePrepare.Start,1),1) = Seq.pulsePrepare.Duration(:,1);
    TX.Amplitude(1:size(Seq.pulsePrepare.Start,1),1) = Seq.pulsePrepare.Amplitude(:,1);
    TX.Frequency(1:size(Seq.pulsePrepare.Start,1),1) = Seq.pulsePrepare.Frequency(:,1);
    TX.Phase(1:size(Seq.pulsePrepare.Start,1),1) = Seq.pulsePrepare.Phase(:,1);
  end

  % excitation pulse
  TX.Start(size(Seq.pulseRefocus.Start,1)+(1:size(Seq.pulseExcitation.Start,1)),(Seq.nTau1>0)+2) = Seq.pulseExcitation.Start - Seq.tEcho/2;
  TX.Duration(size(Seq.pulseRefocus.Start,1)+(1:size(Seq.pulseExcitation.Start,1)),(Seq.nTau1>0)+2) = Seq.pulseExcitation.Duration;
  TX.Amplitude(size(Seq.pulseRefocus.Start,1)+(1:size(Seq.pulseExcitation.Start,1)),(Seq.nTau1>0)+2) = Seq.pulseExcitation.Amplitude;
  TX.Frequency(size(Seq.pulseRefocus.Start,1)+(1:size(Seq.pulseExcitation.Start,1)),(Seq.nTau1>0)+2) = Seq.pulseExcitation.Frequency;
  TX.Phase(size(Seq.pulseRefocus.Start,1)+(1:size(Seq.pulseExcitation.Start,1)),(Seq.nTau1>0)+2) = Seq.pulseExcitation.Phase;

  % refocusing pulse
  TX.Start(1:size(Seq.pulseRefocus.Start,1),(Seq.nTau1>0)+1+(1:Seq.nEcho)) = repmat(Seq.pulseRefocus.Start, 1, Seq.nEcho);
  TX.Duration(1:size(Seq.pulseRefocus.Start,1),(Seq.nTau1>0)+1+(1:Seq.nEcho)) = repmat(Seq.pulseRefocus.Duration, 1, Seq.nEcho);
  TX.Amplitude(1:size(Seq.pulseRefocus.Start,1),(Seq.nTau1>0)+1+(1:Seq.nEcho)) = repmat(Seq.pulseRefocus.Amplitude, 1, Seq.nEcho);
  TX.Frequency(1:size(Seq.pulseRefocus.Start,1),(Seq.nTau1>0)+1+(1:Seq.nEcho)) = repmat(Seq.pulseRefocus.Frequency, 1, Seq.nEcho);
  TX.Phase(1:size(Seq.pulseRefocus.Start,1),(Seq.nTau1>0)+1+(1:Seq.nEcho)) = repmat(Seq.pulseRefocus.Phase, 1, Seq.nEcho) + Seq.refocusingPhaseOffset;

  % sort tRep with excitation pulse
  [TX.Start(:,(Seq.nTau1>0)+2), idx] = sort(TX.Start(:,(Seq.nTau1>0)+2));
  TX.Duration(:,(Seq.nTau1>0)+2) = TX.Duration(idx,(Seq.nTau1>0)+2);
  TX.Amplitude(:,(Seq.nTau1>0)+2) = TX.Amplitude(idx,(Seq.nTau1>0)+2);
  TX.Frequency(:,(Seq.nTau1>0)+2) = TX.Frequency(idx,(Seq.nTau1>0)+2);
  TX.Phase(:,(Seq.nTau1>0)+2) = TX.Phase(idx,(Seq.nTau1>0)+2);

  if Seq.tPreheatRfAmp > 0
    % pre-heat amplifier
    TX.Start(:,1) = [TX.Start(1,1)-Seq.tPreheatRfAmp-50e-6; TX.Start(1:end-1,1)];
    TX.Duration(:,1) = [Seq.tPreheatRfAmp; TX.Duration(1:end-1,1)];
    TX.Amplitude(:,1) = [0; TX.Amplitude(1:end-1,1)];
    TX.Frequency(:,1) = [Seq.fLarmor; TX.Frequency(1:end-1,1)];
    TX.Phase(:,1) = [0; TX.Phase(1:end-1,1)];
  end

  % several trains
  TX.Start = repmat(TX.Start, 1, 1, max(1, Seq.nTau1)+Seq.nTau1SteadyState);
  TX.Duration = repmat(TX.Duration, 1, 1, max(1, Seq.nTau1)+Seq.nTau1SteadyState);
  TX.Amplitude = repmat(TX.Amplitude, 1, 1, max(1, Seq.nTau1)+Seq.nTau1SteadyState);
  TX.Frequency = repmat(TX.Frequency, 1, 1, max(1, Seq.nTau1)+Seq.nTau1SteadyState);
  TX.Phase = repmat(TX.Phase, 1, 1, max(1, Seq.nTau1)+Seq.nTau1SteadyState);

  % adjust data of preparation pulses

  if ~isemptyfield(Seq.pulsePrepare, 'Start') && size(Seq.pulsePrepare.Start, 2) > 1
    TX.Start(1:size(Seq.pulsePrepare.Start,1),1,:) = Seq.pulsePrepare.Start;
  end
  if ~isemptyfield(Seq.pulsePrepare, 'Duration') && size(Seq.pulsePrepare.Duration, 2) > 1
    TX.Duration(1:size(Seq.pulsePrepare.Duration,1),1,:) = Seq.pulsePrepare.Duration;
  end
  if ~isemptyfield(Seq.pulsePrepare, 'Amplitude') && size(Seq.pulsePrepare.Amplitude, 2) > 1
    TX.Amplitude(1:size(Seq.pulsePrepare.Amplitude,1),1,:) = Seq.pulsePrepare.Amplitude;
  end
  if ~isemptyfield(Seq.pulsePrepare, 'Frequency') && size(Seq.pulsePrepare.Frequency, 2) > 1
    TX.Frequency(1:size(Seq.pulsePrepare.Frequency,1),1,:) = Seq.pulsePrepare.Frequency;
  end
  if ~isemptyfield(Seq.pulsePrepare, 'Phase') && size(Seq.pulsePrepare.Phase, 2) > 1
    TX.Phase(1:size(Seq.pulsePrepare.Phase,1),1,:) = Seq.pulsePrepare.Phase;
  end

  TX.Repeat = [];

  numtRepRelax = Seq.nEcho+2+(Seq.nTau1>0);

  %% Damp Coil
  if Seq.DampCoil.Enable
    if HW.TX(iDevice).DampCoil.Enable
      warning('PD:sequence_RecoveryCPMG:DoubleDamp', ...
        'Damping is activated in HW.TX.DampCoil. Are you sure you want to include additional damping on the sequence level?');
    end
    if ~isemptyfield(Seq, 'DigitalIO')
      warning('PD:sequence_RecoveryCPMG:DampCoil', ...
        'Damping Coil not supported if Seq.DigitalIO is already defined')
      Seq.DampCoil.Enable = false;
    else
      % switch Digital IO on after every pulse
      DigitalIODamp.SetTime = nan(2, Seq.nEcho+2+(Seq.nTau1>0));
      DigitalIODamp.SetValue = DigitalIODamp.SetTime;
      DigitalIODamp.SetTime(:,(Seq.nTau1>0)+1+(1:Seq.nEcho)) = ...
        repelem(TX.Start(size(Seq.pulseRefocus.Start,1),(Seq.nTau1>0)+1+(1:Seq.nEcho)) + ...
                TX.Duration(size(Seq.pulseRefocus.Start,1),(Seq.nTau1>0)+1+(1:Seq.nEcho)), 2, 1) + ...
        repmat([0; Seq.DampCoil.Duration], [1,Seq.nEcho]) + Seq.DampCoil.Delay;
      if Seq.DampCoil.DampExcitation
        % Set special duration after 90 degrees excitation pulse
        DigitalIODamp.SetTime(:,(Seq.nTau1>0)+1) = ...
          TX.Start(size(Seq.pulseExcitation.Start,1),(Seq.nTau1>0)+1) + ...
          TX.Duration(size(Seq.pulseExcitation.Start,1),(Seq.nTau1>0)+1) + ...
          [0;Seq.DampCoil.DurationExcitation] + Seq.DampCoil.DelayExcitation;
      end
      DigitalIODamp.SetValue = repmat([2.^(Seq.DampCoil.DOChannel-1); 0], [1,Seq.nEcho+2+(Seq.nTau1>0)]);
      DigitalIODamp.SetValue(isnan(DigitalIODamp.SetTime)) = NaN;
      DigitalIODamp.Repeat = [zeros(1, (Seq.nTau1>0) + 1 + double(Seq.nEcho>1)), ones(1, max(0,Seq.nEcho-1)), 0];
      % several trains
      DigitalIODamp.SetTime = repmat(DigitalIODamp.SetTime, 1, 1, max(1, Seq.nTau1)+Seq.nTau1SteadyState);
      DigitalIODamp.SetValue = repmat(DigitalIODamp.SetValue, 1, 1, max(1, Seq.nTau1)+Seq.nTau1SteadyState);
      DigitalIODamp.Repeat = repmat(DigitalIODamp.Repeat, 1, 1, max(1, Seq.nTau1)+Seq.nTau1SteadyState);
      Seq.DigitalIO = DigitalIODamp;
      Seq.DigitalIO.Repeat = [];
    end
  end


  %% AQ
  Seq.tAQEchoDelay = round(Seq.tAQEchoDelay.*HW.RX(iDevice).fSample)/HW.RX(iDevice).fSample;
  AQ.fSample = Seq.fSample;
  AQ.nSamples = round(Seq.tAQEcho*AQ.fSample);
  if AQ.nSamples==0, AQ.nSamples = 1; end
  Seq.nSampleAQEcho = AQ.nSamples;
  % single CPMG Echo train [preparation, excitation, nEchoes, wait]
  AQ.Start = [NaN(1,Seq.nTau1>0), NaN(1, 1+Seq.SteadyState_PreShots180), ...
    repmat(Seq.tEcho/2-AQ.nSamples/AQ.fSample/2+Seq.tAQEchoDelay, [1, Seq.nEcho-Seq.SteadyState_PreShots180]), NaN];
  if Seq.tAQFID >= 0
    AQ.Start(2,:) = NaN;
    AQ.Start(2,(Seq.nTau1>0)+2) = AQ.Start(1,(Seq.nTau1>0)+2);

    if isempty(Seq.pulseExcitation.Start)
      AQ.Start(1,(Seq.nTau1>0)+2) = Seq.tAQFIDStart - Seq.tRep((Seq.nTau1>0)+1);
    else
      AQ.Start(1,(Seq.nTau1>0)+2) = ...
        max(Seq.tAQFIDStart, ...
            Seq.pulseExcitation.Start(end,min(end,(Seq.nTau1>0)+1)) + ...
            Seq.pulseExcitation.Duration(end,min(end,(Seq.nTau1>0)+1)) + ...
            get_DeadTimeTX2RX(HW, Seq.fSampleFID)) ...
        - Seq.tRep((Seq.nTau1>0)+1);
    end
  end
  if isscalar(Seq.AQPhaseOffset)
    Seq.AQPhaseOffset = repmat(Seq.AQPhaseOffset, ...
      [size(AQ.Start, 1), size(Seq.tRep, 2), size(Seq.tRep, 3)]);
    % AQ.Repeat = repmat([zeros(1,Seq.nTau1>0),0,zeros(1,double(Seq.nEcho>1)),ones(1,max(0,Seq.nEcho-1)),0], 1, max(1,Seq.nTau1)+Seq.nTau1SteadyState);
  end
  if ischar(Seq.AQPhaseOffset)
    if strcmp(Seq.AQPhaseOffset, 'same')
      Seq.AQPhaseOffset = Seq.TXPhaseOffset;
    end
    if strcmp(Seq.AQPhaseOffset, 'alternate')
      Seq.AQPhaseOffset = repmat([zeros(1,Seq.nTau1>0),0,0+90*(-1).^(1:Seq.nEcho),0], ...
        [1, Seq.nTau1+Seq.nTau1SteadyState]);
    end
  end
  if isinf(Seq.tAQFID)
    Seq.tAQFID = Seq.pulseRefocus.Start(1) - AQ.Start(1,(Seq.nTau1>0)+2) - get_DeadTimeRX2TX(HW, AQ.fSample(1));
    if Seq.tAQFID < 1/AQ.fSample(1)
      error('PD:sequence_RecoveryCPMG:NoAQFID', ...
        'tEcho too short or rf pulses too long to acquire FID (or solid echo).');
    end
  end
  % several CPMG Echo trains
  AQ.Start = cat(3, repmat(AQ.Start, [1, 1, Seq.nTau1SteadyState])*nan, ...
    repmat(AQ.Start, [1, 1, max(1, Seq.nTau1)]));
  % extent the acquisition windows to their backs
  AQ.nSamples = round((Seq.tAQEcho)*AQ.fSample);
  if AQ.nSamples==0, AQ.nSamples = 1; end
  Seq.nSampleAQEcho = AQ.nSamples;
  if Seq.ResetDDSPhaseAtExcitation
    % Reset phase of transmission and acquisition at excitation pulse
    % (i.e. after the preparation pulse). This makes sure that the absolute
    % phase of the excitation and inversion pulses are the same for the phase
    % cycle and averaging (even if the preparation time increases).
    AQ.ResetPhases = repmat([0, 1, zeros(1,(Seq.nTau1>0)+Seq.nEcho)], 1, 1, max(1, Seq.nTau1)+Seq.nTau1SteadyState);

    % FIXME: Do we care about the absolute phase of the preparation pulse(s)?
    % AQ.ResetPhases(1,1,:) = 1;
  else
    AQ.ResetPhases = repmat([1, 0, zeros(1,(Seq.nTau1>0)+Seq.nEcho)], 1, 1, max(1, Seq.nTau1)+Seq.nTau1SteadyState);
  end
  AQ.Phase    = Seq.AQPhaseOffset;
  AQ.fSample  = AQ.fSample + zeros(size(AQ.Start));
  AQ.Frequency= AQ.Frequency + zeros(size(AQ.Start));
  AQ.nSamples = AQ.nSamples + zeros(size(AQ.Start));
  if Seq.tAQFID >= 0
    AQ.fSample(1,(Seq.nTau1>0)+2,:) = Seq.fSampleFID;
    Seq.nSampleAQFID = round((Seq.tAQFID)*AQ.fSample(1,(Seq.nTau1>0)+2,:));
    if Seq.nSampleAQFID==0, Seq.nSampleAQFID = 1; end
    AQ.nSamples(1,(Seq.nTau1>0)+2,:) = Seq.nSampleAQFID;
  else
    Seq.nSampleAQFID = 0;
  end


  %% Command load time
  if isempty(Seq.CLTime)
    Seq.CLTime = repmat([3e-6*ones(1,Seq.nTau1>0),3e-6,1.448e-6,0400e-9+zeros(1,Seq.nEcho)], 1, 1, max(1, Seq.nTau1)+Seq.nTau1SteadyState);
  elseif isscalar(Seq.CLTime)
    % use uniform command load time at all tReps
    Seq.CLTime = Seq.CLTime + zeros(size(AQ.Start));
  end

  % Set minimum tOffset of tRep that originally contained the excitation pulse
  % so that it can be "consumed" by tOffset of the subsequent tRep (that now
  % contains the excitation pulse).
  if isemptyfield(TX, 'BlankOffset'), TX.BlankOffset = HW.TX(iDevice).BlankOffset; end
  if isempty(Seq.pulseExcitation.Start)
    % excitation pulse was removed
    % Allow a negative tOffset so that this tRep can be "consumed" by the
    % *preceding* one (e.g., for spin-lock pulse).
    % FIXME: What would be the "best" compromise to allow the preceding and the
    % subsequent tRep to "overlap" into this tRep?
    pulseExcitationStart = -Seq.tEcho/2 + 8e-6;  % Why 8 us???
    if Seq.nTau1 > 0
      pulsePreparationEnd = max(Seq.pulsePrepare.Start(:,1)+Seq.pulsePrepare.Duration(:,1), [], 'omitnan');
      pulseExcitationStart = pulseExcitationStart + (pulsePreparationEnd - Seq.tRep(1));
    end
  else
    pulseExcitationStart = Seq.pulseExcitation.Start(1) - Seq.tEcho/2;
  end
  Seq.tOffset(:,(Seq.nTau1>0)+1,:) = max(Seq.tOffset(:,(Seq.nTau1>0)+1,:), ...
    -Seq.tRep(:,(Seq.nTau1>0)+1,:) - pulseExcitationStart + ...
    sum(Seq.CLTime(:,(Seq.nTau1>0)+1:2,:),2) + TX.BlankOffset + 100/HW.MMRT(iDevice).fSystem);


  %% Grad
  if ~isinf(Seq.thicknessSlice)
    if Seq.tAQFID >= 0
      warning('PD:sequence_RecoveryCPMG:SliceAndFID', ...
        ['Acquiring an FID whilst applying a slice gradient is difficult. ' ...
        'Only "off-slice" signal might be accesible.']);
    end

    gradAmp = Seq.tRefocusBW/(HW.GammaDef/(2*pi) * Seq.thicknessSlice);
    maxAmp = HW.Grad(iDevice).MaxAmp(HW.Grad(iDevice).Channel2xyzB(HW.Grad(iDevice).Slice.channel));
    if gradAmp > maxAmp
      error('PD:sequence_RecoveryCPMG:sliceAmplitudeTooHigh', ...
        'Slice gradient amplitude (%.3f mT/m) too high (max. %.3f mT/m).', ...
        gradAmp*1e3, maxAmp*1e3);
    end

    if Seq.PowerSupply.EnableAnalogChannel > 0
      % set function that configures the external power supply with USB commands

      % keep previous Seq.PreStartPPGFcn
      if isemptyfield(Seq, 'PreStartPPGFcn')
        oldPreStartPPGFcn = [];
      else
        oldPreStartPPGFcn = Seq.PreStartPPGFcn;
      end

      Seq.PreStartPPGFcn = @(sq,nn) SetPowerSupplyAmp(sq,nn, oldPreStartPPGFcn, HW, gradAmp);
    end

    % model for continuous current at constant driving voltage
    % % Calculate the amplitude that would be reached after infinite time with
    % % constant driving voltage where the desired amplitude is reached after
    % % HW.Grad(iDevice).Slice.tRamp.
    % % IEff = Igrad / (1 - exp(-R/L*tRamp))
    % ampEff = gradAmp / ...
    %   (1 - exp(-HW.Grad(iDevice).LoadRin(HW.Grad(iDevice).Slice.channel) / ...
    %             HW.Grad(iDevice).Inductance(HW.Grad(iDevice).Slice.channel) * ...
    %             HW.Grad(iDevice).Slice.tRamp));
    %
    % % calculate amplitude for steps
    % C = (1-gradAmp/ampEff)^(1/Seq.Slice.nRamp)/(1-(1-gradAmp/ampEff)^(1/Seq.Slice.nRamp));
    %
    % gradAmpStep = zeros(Seq.Slice.nRamp+1, 1);
    % for n = 1:Seq.Slice.nRamp
    %   gradAmpStep(n+1) = (ampEff + C*gradAmpStep(n))/(1+C);
    % end

    % Model for discretized linear slopes where the voltage at the end of each
    % segment is the same as at the end of the previous segment.
    % V[n] = sum_(m=1)^(n) (delta_I[n]) * R  + L*delta_I[n]/delta_t
    % V[n] = V[n-1]
    % => delta_I[n] = delta_I[1] / (R/L*delta_t + 1)^(n-1)
    % I[n] = I[1] * sum_(m=0)^(n) (1/(R/L*delta_t + 1)^m)
    %   X = 1/(R/L*delta_t + 1):
    % I[n] = I[1] * (1-X^n)/(1-X)
    % => for n=nRamp (I[nRamp] = Igrad)
    % I[1] = Igrad * (1-X)/(1-X^nRamp)
    X = 1 / ...
      (HW.Grad(iDevice).LoadRin(HW.Grad(iDevice).Slice.channel) / ...
       HW.Grad(iDevice).Inductance(HW.Grad(iDevice).Slice.channel) * ...
       HW.Grad(iDevice).Slice.tRamp/Seq.Slice.nRamp + 1);
    amp1 = gradAmp / ( (1 - X^Seq.Slice.nRamp) / (1 - X) );
    gradAmpStep = [0; amp1*(1-X.^(1:Seq.Slice.nRamp).')/(1-X)];

    % FIXME: Do we need a different calculation for a "best" shape for ramping
    % down?

    % Check if time between CPMGs is long enough for Seq.Slice.offBetweenCPMG
    minIETime = 1.5*HW.Grad(iDevice).Slice.tRamp + HW.Grad(iDevice).Slice.tEC + Seq.CLTime(numtRepRelax:numtRepRelax:end);
    if any(Seq.tRep(numtRepRelax:numtRepRelax:end) < minIETime)
      error('PD:sequence_Recovery_CPMG:RepetitionRateTooShort', ...
        'The time between echo trains is too short to turn off the gradients. It must be longer than %.2f%sms.', ...
        1e3*min(minIETime), 160);
      % Seq.Slice.offBetweenCPMG = false;
    end

    GradSlice.Time = nan(1 + Seq.Slice.nRamp + 2*Seq.Slice.offAfterPreparation, size(AQ.Start, 2), size(AQ.Start, 3));
    GradSlice.Repeat = ones(1, size(AQ.Start, 2), size(AQ.Start, 3));
    GradSlice.Amp = GradSlice.Time;

    if Seq.PowerSupply.EnableAnalogChannel > 0
      DigiPS.SetTime = nan(2, size(AQ.Start, 2), size(AQ.Start, 3));
      DigiPS.SetValue = DigiPS.SetTime;
    end

    % ramp up
    if Seq.Slice.offBetweenCPMG
      GradSlice.Time(1:(1+Seq.Slice.nRamp),1,:) = repmat(linspace(-HW.Grad(iDevice).Slice.tRamp, 0, Seq.Slice.nRamp+1), 1, 1, size(GradSlice.Time, 3)) - HW.Grad(iDevice).Slice.tEC;
      if Seq.nTau1 > 0
        GradSlice.Time(1:(1+Seq.Slice.nRamp),1,:) = GradSlice.Time(1:(1+Seq.Slice.nRamp),1,:) + Seq.pulsePrepare.Start(1);
      else
        GradSlice.Time(1:(1+Seq.Slice.nRamp),1,:) = GradSlice.Time(1:(1+Seq.Slice.nRamp),1,:) + Seq.pulseExcitation.Start(1);
      end
      GradSlice.Amp(1:(1+Seq.Slice.nRamp),1,:) = repmat(gradAmpStep, 1, 1, size(GradSlice.Amp, 3));
      GradSlice.Repeat(1,1,:) = 0;
      if Seq.PowerSupply.EnableAnalogChannel > 0
        DigiPS.SetTime(1,1,:) = GradSlice.Time(1,1,:) - Seq.PowerSupply.EnableAnalogTOffset;
        DigiPS.SetValue(1,1,:) = 2^(Seq.PowerSupply.EnableAnalogChannel-1);
      end
    else
      GradSlice.Time(1:(1+Seq.Slice.nRamp),1) = linspace(-HW.Grad(iDevice).Slice.tRamp, 0, Seq.Slice.nRamp+1) - HW.Grad(iDevice).Slice.tEC;
      if Seq.nTau1 > 0
        GradSlice.Time(1:(1+Seq.Slice.nRamp),1) = GradSlice.Time(1:(1+Seq.Slice.nRamp),1) + Seq.pulsePrepare.Start(1);
      else
        GradSlice.Time(1:(1+Seq.Slice.nRamp),1) = GradSlice.Time(1:(1+Seq.Slice.nRamp),1) + Seq.pulseExcitation.Start(1);
      end
      GradSlice.Amp(1:(1+Seq.Slice.nRamp),1) = gradAmpStep;
      GradSlice.Repeat(1) = 0;
      if Seq.PowerSupply.EnableAnalogChannel > 0
        DigiPS.SetTime(1,1) = GradSlice.Time(1,1) - Seq.PowerSupply.EnableAnalogTOffset;
        DigiPS.SetValue(1,1) = 2^(Seq.PowerSupply.EnableAnalogChannel-1);
      end
    end

    if Seq.nTau1 > 0 && Seq.Slice.offAfterPreparation
      % Check if gap between preparation and excitation is long enough
      selectedTau1 = Seq.tRep(1,1,:) - (Seq.pulsePrepare.Start(1) + sum(Seq.pulsePrepare.Duration)) + Seq.pulseExcitation.Start(1) ...
        > 2*HW.Grad(iDevice).Slice.tRamp + HW.Grad(iDevice).Slice.tEC + Seq.CLTime(1,1,:) + HW.Grad(iDevice).TimeDelay(HW.Grad(iDevice).Slice.channel);
      % ramp down
      GradSlice.Time((1+Seq.Slice.nRamp)+(1:(1+Seq.Slice.nRamp)),1,selectedTau1) = repmat(linspace(0, HW.Grad(iDevice).Slice.tRamp/2, Seq.Slice.nRamp+1).', 1, 1, sum(selectedTau1)) ...
        + Seq.pulsePrepare.Start(1) + sum(Seq.pulsePrepare.Duration) + HW.Grad(iDevice).TimeDelay(HW.Grad(iDevice).Slice.channel);
      GradSlice.Amp((1+Seq.Slice.nRamp)+(1:(1+Seq.Slice.nRamp)),1,selectedTau1) = repmat(gradAmp-gradAmpStep, 1, 1, sum(selectedTau1));
      % ramp up
      GradSlice.Time(1:(1+Seq.Slice.nRamp),2,selectedTau1) = repmat(linspace(-HW.Grad(iDevice).Slice.tRamp, 0, Seq.Slice.nRamp+1).', 1, 1, sum(selectedTau1)) ...
         - HW.Grad(iDevice).Slice.tEC + Seq.pulseExcitation.Start(1);
      GradSlice.Amp(1:(1+Seq.Slice.nRamp),2,selectedTau1) = repmat(gradAmpStep, 1, 1, sum(selectedTau1));
      GradSlice.Repeat(1,1:(1+Seq.Slice.nRamp),selectedTau1) = 0;
      if Seq.PowerSupply.EnableAnalogChannel > 0
        selectedTau1PSPrep = Seq.tRep(1,1,:) - (Seq.pulsePrepare.Start(1) + sum(Seq.pulsePrepare.Duration)) + Seq.pulseExcitation.Start(1) ...
          > 2*HW.Grad(iDevice).Slice.tRamp + HW.Grad(iDevice).Slice.tEC + Seq.CLTime(1,1,:) + HW.Grad(iDevice).TimeDelay(HW.Grad(iDevice).Slice.channel) ...
          + Seq.PowerSupply.EnableAnalogTPostset + Seq.PowerSupply.EnableAnalogTOffset;
        DigiPS.SetTime(2,1,selectedTau1PSPrep) = GradSlice.Time(2*(1+Seq.Slice.nRamp),1) + Seq.PowerSupply.EnableAnalogTPostset;
        DigiPS.SetValue(2,1,selectedTau1PSPrep) = 0;
        DigiPS.SetTime(1,2,selectedTau1PSPrep) = GradSlice.Time(1,1) - Seq.PowerSupply.EnableAnalogTOffset;
        DigiPS.SetValue(1,2,selectedTau1PSPrep) = 2^(Seq.PowerSupply.EnableAnalogChannel-1);
      end
    else
      selectedTau1 = false(1, 1, size(Seq.tRep, 3));
    end

    % keep gradient amp during tReps
    if Seq.Slice.offBetweenCPMG
      GradSlice.Time(1,3:end-1,:) = 0;
      GradSlice.Amp(1,3:end-1,:) = gradAmp;
      GradSlice.Time(1,2,~selectedTau1) = 0;
      GradSlice.Amp(1,2,~selectedTau1) = gradAmp;
    else
      % FIXME: Do not overwrite the tReps with Seq.Slice.offAfterPreparation
      GradSlice.Time(1,2:end-1) = 0;
      GradSlice.Amp(1,2:end-1) = gradAmp;
    end

    % ramp down
    if Seq.Slice.offBetweenCPMG
      GradSlice.Time(1:(1+Seq.Slice.nRamp),end,:) = repmat(linspace(0, HW.Grad(iDevice).Slice.tRamp/2, Seq.Slice.nRamp+1), ...
        1, 1, size(GradSlice.Time, 3)) + HW.Grad(iDevice).TimeDelay(HW.Grad(iDevice).Slice.channel);
      GradSlice.Amp(1:(1+Seq.Slice.nRamp),end,:) = repmat(gradAmp-gradAmpStep, 1, 1, size(GradSlice.Amp, 3));
      GradSlice.Repeat(1,end,:) = 0;
      if Seq.PowerSupply.EnableAnalogChannel > 0
        DigiPS.SetTime(1,end,:) = GradSlice.Time(1+Seq.Slice.nRamp,end,:) + Seq.PowerSupply.EnableAnalogTPostset;
        DigiPS.SetValue(1,end,:) = 0;
        crossingTau1PS = Seq.tRep(1,end,:) + DigiPS.SetTime(1,1,:) ...
          <= DigiPS.SetTime(1,end,:);
      end
    else
      GradSlice.Time(1:(1+Seq.Slice.nRamp),end) = linspace(0, HW.Grad(iDevice).Slice.tRamp/2, Seq.Slice.nRamp+1) + ...
        HW.Grad(iDevice).TimeDelay(HW.Grad(iDevice).Slice.channel);
      GradSlice.Amp(1:(1+Seq.Slice.nRamp),end) = gradAmp-gradAmpStep;
      GradSlice.Repeat(end) = 0;
      if Seq.PowerSupply.EnableAnalogChannel > 0
        DigiPS.SetTime(1,end,end) = GradSlice.Time(1+Seq.Slice.nRamp,end) + Seq.PowerSupply.EnableAnalogTPostset;
        DigiPS.SetValue(1,end,end) = 0;
        crossingTau1PS = false(1,1,size(DigiPS.SetTime,3));
      end
    end

    % FIXME: This might overwrite previously set gradients. Might be a non-issue
    % if the slice gradient is on a different channel than the imaging gradients.
    if isscalar(Grad(HW.Grad(iDevice).Channel2xyzB(HW.Grad(iDevice).Slice.channel)).Time) && ...
        isnan(Grad(HW.Grad(iDevice).Channel2xyzB(HW.Grad(iDevice).Slice.channel)).Time)
      Grad(HW.Grad(iDevice).Channel2xyzB(HW.Grad(iDevice).Slice.channel)).Time = GradSlice.Time;
      Grad(HW.Grad(iDevice).Channel2xyzB(HW.Grad(iDevice).Slice.channel)).Amp = GradSlice.Amp;
    else
      Grad(HW.Grad(iDevice).Channel2xyzB(HW.Grad(iDevice).Slice.channel)).Time(~isnan(GradSlice.Time)) = GradSlice.Time(~isnan(GradSlice.Time));
      Grad(HW.Grad(iDevice).Channel2xyzB(HW.Grad(iDevice).Slice.channel)).Amp(~isnan(GradSlice.Amp)) = GradSlice.Amp(~isnan(GradSlice.Amp));
    end
    Grad(HW.Grad(iDevice).Channel2xyzB(HW.Grad(iDevice).Slice.channel)).Repeat(GradSlice.Repeat==0) = 0;
    Grad(HW.Grad(iDevice).Channel2xyzB(HW.Grad(iDevice).Slice.channel)).Repeat = [];
  end


  if Seq.ResetDDSPhaseAtExcitation && ~isemptyfield(Seq.pulseExcitation, 'tCenterOffset')
    % If the pulse shape function returned a value for the center of the rf
    % pulses (e.g., Pulse_Rect_SolidEcho), move the phase reference of the DDS
    % to that time.

    % FIXME: Would we like to move the reference point for the DDS phase also
    %        for pulses that don't return "tCenterOffset" (e.g., by assuming 0)?

    % FIXME: The following is probably not correct if we'd start to support
    %        off-center (PRESS) inversion pulses.
    TX.Phase(:,2:end,:) = mod(bsxfun(@minus, ...
      TX.Phase(:,2:end,:), ...
      -TX.Frequency(1,2,:) * Seq.pulseExcitation.tCenterOffset * 360), 360);
    AQ.Phase(:,2:end,:) = mod(bsxfun(@minus, ...
      AQ.Phase(:,2:end,:), ...
      -AQ.Frequency(1,2,:) * Seq.pulseExcitation.tCenterOffset * 360), 360);
  end


  %% User function hook for sequence manipulation
  Seq.T2Prepare.tRepsPrepare = (Seq.nTau1>0);
  if ~isempty(Seq.Function_Prepare_Measurement)
    [HW, Seq, AQ, TX, Grad] = Seq.Function_Prepare_Measurement(HW, Seq, AQ, TX, Grad);
  end


  %% Averages (phase cycle)
  if Seq.nTau1>0
    % from increment to accumulated phase offset
    tempTXPhasePrepareIncrement = cumsum(repmat(permute(Seq.SeqAverage.TXPhasePrepareIncrement, [1,4,2,3]), ...
      [size(Seq.pulsePrepare.Start,1), Seq.T2Prepare.tRepsPrepare, Seq.SeqAverage.average/numel(Seq.SeqAverage.TXPhasePrepareIncrement), max(1, Seq.nTau1)]), 3) - Seq.SeqAverage.TXPhasePrepareIncrement(1);
    % apply to all CPMG trains
    tempTXPhasePrepareOffset = repmat(permute(Seq.SeqAverage.TXPhasePrepareOffset, [1,4,2,3]), ...
      [size(Seq.pulsePrepare.Start,1), Seq.T2Prepare.tRepsPrepare, 1, max(1, Seq.nTau1)]);
  else
    tempTXPhasePrepareIncrement = [];
    tempTXPhasePrepareOffset = [];
  end
  % from increment to accumulated phase offset
  tempTXPhaseExcitationIncrement = cumsum(repmat(permute(Seq.SeqAverage.TXPhaseExcitationIncrement, [1,4,2,3]), ...
    [size(Seq.pulseExcitation.Start,1), 1, Seq.SeqAverage.average/numel(Seq.SeqAverage.TXPhaseExcitationIncrement), max(1, Seq.nTau1)]), 3) - Seq.SeqAverage.TXPhaseExcitationIncrement(1);
  tempTXPhaseRefocusIncrement = cumsum(repmat(permute(Seq.SeqAverage.TXPhaseRefocusIncrement, [1,4,2,3]), ...
    [size(Seq.pulseRefocus.Start,1), Seq.nEcho+1, Seq.SeqAverage.average/numel(Seq.SeqAverage.TXPhaseRefocusIncrement), max(1, Seq.nTau1)]), 3) - Seq.SeqAverage.TXPhaseRefocusIncrement(1);
  tempAQPhaseIncrement = cumsum(repmat(permute(Seq.SeqAverage.AQPhaseIncrement, [1,4,2,3]), ...
    [size(AQ.Start, 1), 2+Seq.T2Prepare.tRepsPrepare+Seq.nEcho, Seq.SeqAverage.average/numel(Seq.SeqAverage.AQPhaseIncrement), max(1, Seq.nTau1)]), 3) - Seq.SeqAverage.AQPhaseIncrement(1);

  % apply offset to all CPMG trains
  tempTXPhaseExcitationOffset = repmat(permute(Seq.SeqAverage.TXPhaseExcitationOffset, [1,4,2,3]), ...
    [size(Seq.pulseExcitation.Start,1), 1, 1, max(1, Seq.nTau1)]);
  tempTXPhaseRefocusOffset = repmat(permute(Seq.SeqAverage.TXPhaseRefocusOffset, [1,4,2,3]), ...
    [size(Seq.pulseRefocus.Start,1), Seq.nEcho+1, 1, max(1, Seq.nTau1)]);
  tempAQPhaseOffset = repmat(permute(Seq.SeqAverage.TXPhaseExcitationOffset, [1,4,2,3]), ...
    [size(AQ.Start, 1), 2+Seq.T2Prepare.tRepsPrepare+Seq.nEcho, 1, max(1, Seq.nTau1)]);

  numelPulsePrepare = size(tempTXPhasePrepareIncrement, 1);
  numelPulseExcite = size(tempTXPhaseExcitationIncrement, 1);
  numelPulseRefocus = size(tempTXPhaseRefocusIncrement, 1);
  maxNumelPulse = max([numelPulsePrepare, numelPulseExcite, numelPulseRefocus]);
  tempTXPhaseAllOffset = cat(2, ...
    cat(1, tempTXPhasePrepareIncrement+tempTXPhasePrepareOffset, zeros(maxNumelPulse-numelPulsePrepare, size(tempTXPhasePrepareIncrement, 2), Seq.SeqAverage.average, max(1, Seq.nTau1))), ...
    cat(1, tempTXPhaseExcitationIncrement+tempTXPhaseExcitationOffset, zeros(maxNumelPulse-numelPulseExcite, size(tempTXPhaseExcitationIncrement, 2), Seq.SeqAverage.average, max(1, Seq.nTau1))), ...
    cat(1, tempTXPhaseRefocusIncrement+tempTXPhaseRefocusOffset, zeros(maxNumelPulse-numelPulseRefocus, size(tempTXPhaseRefocusIncrement, 2), Seq.SeqAverage.average, max(1, Seq.nTau1))));
  tempTXPhaseAllOffset = reshape(tempTXPhaseAllOffset, maxNumelPulse, Seq.nEcho+2+Seq.T2Prepare.tRepsPrepare, []);
  % excitation pulse is first pulse in the same tRep as the first refocus pulse
  tempTXPhaseAllOffset(numelPulseExcite+(1:numelPulseRefocus),size(tempTXPhasePrepareIncrement,2)+2,:) = ...
    tempTXPhaseAllOffset(1:numelPulseRefocus,size(tempTXPhasePrepareIncrement,2)+2,:);
  tempTXPhaseAllOffset(1:numelPulseExcite,size(tempTXPhasePrepareIncrement,2)+2,:) = ...
    tempTXPhaseAllOffset(1:numelPulseExcite,size(tempTXPhasePrepareIncrement,2)+1,:);
  % "phase cycle" for excitation pulse (e.g. for solid echo pulse)
  PulseExcitation.Phase = Seq.SeqAverage.TXPhaseExcitationIncrement;
  % only applies to solid echo pulses
  PulseExcitation.PhaseSolidInversionOffset = -90;
  pulseExcitationPhC = Seq.excitationPulse(HW, 0, PulseExcitation);
  tempTXPhaseAllOffset(1:numelPulseExcite,size(tempTXPhasePrepareIncrement,2)+2,2:2:end) = ...
    repmat(pulseExcitationPhC.Phase-Seq.pulseExcitation.Phase, [1, 1, floor(size(tempTXPhaseAllOffset, 3)/2)]);
  tempTXPhaseAllOffset(1:numelPulseExcite,size(tempTXPhasePrepareIncrement,2)+1,:) = 0;
  tempAQPhaseAllOffset = reshape(tempAQPhaseIncrement+tempAQPhaseOffset, ...
    size(AQ.Start, 1), Seq.nEcho+2+Seq.T2Prepare.tRepsPrepare, []);

  % repeat pulse programs (and apply phase offset) - but not the pre-shots
  Seq.tRep = reshape(cat(3, [Seq.tRep(:,1:end-1,1:Seq.nTau1SteadyState), Seq.tRep(:,end,1:Seq.nTau1SteadyState)+Seq.SeqAverage.averageBreak], ... % pre-shots
    reshape(permute(repmat([Seq.tRep(:,1:end-1,(Seq.nTau1SteadyState+1):end), Seq.tRep(:,end,(Seq.nTau1SteadyState+1):end)+Seq.SeqAverage.averageBreak], [1, 1, 1, Seq.SeqAverage.average]), [1 2 4 3]), ...
    size(Seq.tRep,1), size(Seq.tRep,2), [])), ...
    1, []);
  Seq.tOffset = reshape(cat(3, [Seq.tOffset(:,1:end-1,1:Seq.nTau1SteadyState), Seq.tOffset(:,end,1:Seq.nTau1SteadyState)+Seq.SeqAverage.averageBreak], ... % pre-shots
    reshape(permute(repmat([Seq.tOffset(:,1:end-1,(Seq.nTau1SteadyState+1):end), Seq.tOffset(:,end,(Seq.nTau1SteadyState+1):end)+Seq.SeqAverage.averageBreak], [1, 1, 1, Seq.SeqAverage.average]), [1 2 4 3]), ...
    size(Seq.tOffset,1), size(Seq.tOffset,2), [])), ...
    1, []);
  TX.Start = reshape(cat(3, TX.Start(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
    reshape(repmat(TX.Start(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(TX.Start, 1), size(TX.Start, 2), [])), ...
    size(TX.Start, 1), []);
  TX.Duration = reshape(cat(3, TX.Duration(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
    reshape(repmat(TX.Duration(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(TX.Duration, 1), size(TX.Duration, 2), [])), ...
    size(TX.Start, 1), []);
  TX.Frequency = reshape(cat(3, TX.Frequency(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
    reshape(repmat(TX.Frequency(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(TX.Frequency, 1), size(TX.Frequency, 2), [])), ...
    size(TX.Start, 1), []);
  TX.Amplitude = reshape(cat(3, TX.Amplitude(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
    reshape(repmat(TX.Amplitude(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(TX.Amplitude, 1), size(TX.Amplitude, 2), [])), ...
    size(TX.Start, 1), []);
  TX.Phase = reshape(bsxfun(@plus, cat(3, TX.Phase(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
    reshape(repmat(TX.Phase(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(TX.Phase, 1), size(TX.Phase, 2), [])), ...
    reshape(cat(3, repmat(tempTXPhaseAllOffset(:,:,1), [1, 1, Seq.nTau1SteadyState]), ... % pre-shots
    tempTXPhaseAllOffset), [], 2+1*Seq.T2Prepare.tRepsPrepare+Seq.nEcho, max(1, Seq.nTau1)*Seq.SeqAverage.average + Seq.nTau1SteadyState)), ...
    size(TX.Start, 1), []);
  if ~isemptyfield(TX, 'Repeat')
    TX.Repeat = reshape(cat(3, TX.Repeat(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
      repmat(TX.Repeat(:,:,(Seq.nTau1SteadyState+1):end), [1, 1, Seq.SeqAverage.average])), ...
      size(TX.Start, 1), []);
  end

  % Damp Coil
  if Seq.DampCoil.Enable
    Seq.DigitalIO.SetTime  = reshape(cat(3, Seq.DigitalIO.SetTime(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
      reshape(repmat(Seq.DigitalIO.SetTime(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(Seq.DigitalIO.SetTime, 1), size(Seq.DigitalIO.SetTime, 2), [])), ...
      size(Seq.DigitalIO.SetTime, 1), []);
    Seq.DigitalIO.SetValue = reshape(cat(3, Seq.DigitalIO.SetValue(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
      reshape(repmat(Seq.DigitalIO.SetValue(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(Seq.DigitalIO.SetValue, 1), size(Seq.DigitalIO.SetValue, 2), [])), ...
      size(Seq.DigitalIO.SetTime, 1), []);
    %Seq.DigitalIO.Repeat   = repmat(Seq.DigitalIO.Repeat,   1, Seq.SeqAverage.average);
    Seq.DigitalIO.Repeat = [];
  end

  AQ.Start = reshape(cat(3, AQ.Start(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
    reshape(repmat(AQ.Start(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(AQ.Start, 1), size(AQ.Start, 2), [])), ...
    size(AQ.Start, 1), []);
  AQ.fSample = reshape(cat(3, AQ.fSample(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
    reshape(repmat(AQ.fSample(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(AQ.fSample, 1), size(AQ.fSample, 2), [])), ...
    size(AQ.Start, 1), []);
  AQ.Phase = reshape(cat(3, AQ.Phase(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
    reshape(repmat(AQ.Phase(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(AQ.Phase, 1), size(AQ.Phase, 2), [])) + ...
    reshape(cat(3, repmat(tempAQPhaseAllOffset(:,:,1), [1, 1, Seq.nTau1SteadyState]), ... % pre-shots
                   tempAQPhaseAllOffset), ...
            [], 2+1*Seq.T2Prepare.tRepsPrepare+Seq.nEcho, max(1, Seq.nTau1)*Seq.SeqAverage.average + Seq.nTau1SteadyState), ...
    size(AQ.Start, 1), []);
  AQ.ResetPhases = reshape(cat(3, AQ.ResetPhases(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
    reshape(repmat(AQ.ResetPhases(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(AQ.ResetPhases, 1), size(AQ.ResetPhases, 2), [])), ...
    1, []);
  AQ.nSamples = reshape(cat(3, AQ.nSamples(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
    reshape(repmat(AQ.nSamples(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(AQ.nSamples, 1), size(AQ.nSamples, 2), [])), ...
    size(AQ.Start, 1), []);
  if ~isemptyfield(AQ, 'Repeat')
    AQ.Repeat = reshape(cat(3, AQ.Repeat(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
      reshape(repmat(AQ.Repeat(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(AQ.Repeat, 1), size(AQ.Repeat, 2), [])), ...
      size(AQ.Start, 1), []);
  end
  if ~isemptyfield(Seq, 'CLTime')
    Seq.CLTime = reshape(cat(3, Seq.CLTime(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
      reshape(repmat(Seq.CLTime(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(Seq.CLTime, 1), size(Seq.CLTime, 2), [])), ...
      1, []);
  end

  if size(AQ.Frequency, 2) * size(AQ.Frequency, 3)*Seq.SeqAverage.average < size(AQ.Start, 2)
    % FIXME: Should we check the size (multiples) more exactly here?
    AQ.Frequency = repmat(AQ.Frequency(:,1,:), 1, size(AQ.Start,2)./Seq.SeqAverage.average);
  end
  AQ.Frequency = reshape(cat(3, AQ.Frequency(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
    repmat(AQ.Frequency(:,:,(Seq.nTau1SteadyState+1):end), [1, 1, Seq.SeqAverage.average])), ...
    size(AQ.Start, 1), []);

  for t = 1:HW.Grad(iDevice).n
    if ~isscalar(Grad(t).Time)
      if size(Grad(t).Time, 2) == 1, Grad(t).Time = repmat(Grad(t).Time(:,1), 1, size(Seq.tRep,2)); end
      if size(Grad(t).Amp, 2) == 1,  Grad(t).Amp  = repmat(Grad(t).Amp(:,1),  1, size(Seq.tRep,2)); end

      Grad(t).Time = reshape(cat(3, Grad(t).Time(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
        reshape(repmat(Grad(t).Time(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(Grad(t).Time, 1), size(Grad(t).Time, 2), [])), ...
        size(Grad(t).Time, 1), []);
      Grad(t).Amp  = reshape(cat(3, Grad(t).Amp(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
         reshape(repmat(Grad(t).Amp(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(Grad(t).Amp, 1), size(Grad(t).Amp, 2), [])), ...
        size(Grad(t).Amp, 1), []);
      if ~isemptyfield(Grad(t), 'Repeat')
        Grad(t).Repeat = reshape(cat(3, Grad(t).Repeat(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
          reshape(repmat(Grad(t).Repeat(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(Grad(t).Repeat, 1), size(Grad(t).Repeat, 2), [])), ...
          size(Grad(t).Repeat, 1), []);
      end
      if ~Seq.Slice.offBetweenCPMG
        % don't turn off slice gradient between averages
        Grad(t).Time(:,[numtRepRelax:numtRepRelax:end-1,numtRepRelax+1:numtRepRelax:end-1]) = NaN;
        Grad(t).Time(1,[numtRepRelax:numtRepRelax:end-1,numtRepRelax+1:numtRepRelax:end-1]) = HW.Grad(iDevice).TimeDelay(HW.Grad(iDevice).Slice.channel);
        Grad(t).Amp(1,[numtRepRelax:numtRepRelax:end-1,numtRepRelax+1:numtRepRelax:end-1]) = gradAmp;
        Grad(t).Amp(isnan(Grad(t).Time)) = NaN;
      end
    end
  end

  % external power supply analog remote control
  if Seq.PowerSupply.EnableAnalogChannel > 0 && ~isinf(Seq.thicknessSlice)
    crossingTau1PS2 = false(1,size(DigiPS.SetTime,2),size(crossingTau1PS,3));
    crossingTau1PS2(1,1,:) = crossingTau1PS;
    crossingTau1PS2(1,end,:) = crossingTau1PS;
    DigiPS.SetTime = reshape(cat(3, DigiPS.SetTime(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
      reshape(repmat(DigiPS.SetTime(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(DigiPS.SetTime, 1), size(DigiPS.SetTime, 2), [])), ...
      size(DigiPS.SetTime, 1), []);
    DigiPS.SetValue = reshape(cat(3, DigiPS.SetValue(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
      reshape(repmat(DigiPS.SetValue(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(DigiPS.SetValue, 1), size(DigiPS.SetValue, 2), [])), ...
      size(DigiPS.SetValue, 1), []);
    crossingTau1PS2 = reshape(cat(3, crossingTau1PS2(:,:,1:Seq.nTau1SteadyState), ... % pre-shots
      reshape(repmat(crossingTau1PS2(:,:,(Seq.nTau1SteadyState+1):end), [1, Seq.SeqAverage.average, 1]), size(crossingTau1PS2, 1), size(crossingTau1PS2, 2), [])), ...
      size(crossingTau1PS2, 1), []);
    crossingTau1PS2([1,end]) = false;
    DigiPS.SetTime(:,crossingTau1PS2) = NaN;
    DigiPS.SetValue(:,crossingTau1PS2) = NaN;

    allNaN = all(isnan(DigiPS.SetTime), 2);
    DigiPS.SetTime(allNaN,:) = [];
    DigiPS.SetValue(allNaN,:) = [];
    if isemptyfield(Seq, {'DigitalIO', 'SetTime'})
      Seq.DigitalIO = DigiPS;
    else
      Seq.DigitalIO = add_DigitalIO(Seq.DigitalIO, DigiPS);
    end
  end


  %% Get data directly after last AQ window
  AQ.GetData = zeros(size(Seq.tRep));
  if Seq.GetDataEarly
    % also consider measurements without any echoes (only FID)
    AQ.GetData(end-(Seq.nEcho > 1)) = 1;
  end
  if Seq.SeqAverage.GetDataAfterAverages
    % Mark the last acquisition of the last echo train that belongs to the same
    % preparation time.
    AQ.GetData(tRepsPerCMPG*(Seq.nTau1SteadyState+Seq.SeqAverage.average)-1:tRepsPerCMPG*Seq.SeqAverage.average:end) = 1;
  end


  %% sample lift
  if Seq.Lift.useLift && ~Seq.Lift.skipLift
    DigitalIOLift.SetTime = NaN(2, numel(Seq.tRep));
    DigitalIOLift.SetValue = DigitalIOLift.SetTime;
    DigitalIOLift.SetTime(:,end) = [0; 5e-3;];
    DigitalIOLift.SetValue(:,end) = [2^(Seq.Lift.startSetChannel-1) + ...
      bitget(Seq.Lift.useSet-1, 1) * 2^(Seq.Lift.setBit0Channel-1); 0];
    if ~Seq.SubtractEmptyMeasurement
      if isemptyfield(Seq, 'DigitalIO') || isemptyfield(Seq.DigitalIO, 'SetTime')
        Seq.DigitalIO = DigitalIOLift;
      else
        Seq.DigitalIO = add_DigitalIO(Seq.DigitalIO, DigitalIOLift);
      end
    end
    Seq.DigitalIO.Repeat = [];
  end


  %% shorten last tRep
  if Seq.Lift.useLift && ~Seq.Lift.skipLift
    Seq.tRep(end) = DigitalIOLift.SetTime(end) + Seq.CLTime(end);
  else
    % shorten last tRep
    Seq.tRep(end) = Seq.tEcho;
  end

  if ~isinf(Seq.thicknessSlice)
    if Seq.PowerSupply.EnableAnalogChannel > 0
      Seq.tRep(end) = max(Seq.tRep(end), DigiPS.SetTime(1,end) + Seq.CLTime(end), 'omitnan');
    end
    Seq.tRep(end) = max(Seq.tRep(end), ...
      Grad(HW.Grad(iDevice).Channel2xyzB(HW.Grad(iDevice).Slice.channel)).Time(1+Seq.Slice.nRamp,end) + Seq.CLTime(end));
  end

  Seq.TimeToNextSequence = Seq.tSlice - Seq.tRep(end);
  % FIXME: This should use Seq.tOffset(1). But the offsets aren't known yet.
  %        Use an estimated value instead.
  if isempty(Seq.pulseExcitation.Start)
    excStart = 0;
  else
    excStart = Seq.pulseExcitation.Start(1);
  end
  maybe_tOffset1 = ~isinf(Seq.thicknessSlice)*(HW.Grad(iDevice).Slice.tRamp + HW.Grad(iDevice).Slice.tEC) ...
    - excStart;
  if Seq.TimeToNextSequence < maybe_tOffset1
    % FIXME: What do we have to check if there is a preparation pulse
    warning('PD:sequence_RecoveryCPMG:TimeToNextSequenceTooShort', ...
      ['TimeToNextSequence is too short: %.3f s\n', ...
      'It should be longer than %.3f s plus the time needed for post-processing.'], ...
      Seq.TimeToNextSequence, maybe_tOffset1);
    Seq.TimeToNextSequence = maybe_tOffset1 + 20e-3;
  end

  Seq.averageRepetitionTime = sum(Seq.tRep(1:end-1)) + Seq.tRelax;
  Seq.averageBreak = [];
  if sum(Seq.tRep) + maybe_tOffset1 > Seq.averageRepetitionTime
    Seq.averageRepetitionTime = [];
    Seq.averageBreak = maybe_tOffset1;
  end


  %% prepare the sequence
  StartSequence = Seq.StartSequence;
  PollPPGfast = Seq.PollPPGfast;
  GetRawData = Seq.GetRawData;
  PostProcessSequence = Seq.PostProcessSequence;
  Seq.StartSequence = 0;
  Seq.PollPPGfast = 0;
  Seq.GetRawData = 0;
  Seq.PostProcessSequence = 0;

  [~, Seq] = set_sequence(HW, Seq, AQ, TX, Grad);

  Seq.StartSequence = StartSequence;
  Seq.PollPPGfast = PollPPGfast;
  Seq.GetRawData = GetRawData;
  Seq.PostProcessSequence = PostProcessSequence;
  Seq.PreProcessSequence = 0;


  if Seq.TimeToNextSequence < Seq.tOffset(1) + Seq.CLTime(1)  % FIXME: What is a better minimum value for TimeToNextSequence?
    warning('PD:sequence_RecoveryCPMG:TimeToNextSequenceTooShort', ...
      ['TimeToNextSequence is too short: %.3f s\n', ...
      'It should be longer than %.3f s plus the time needed for post-processing.'], ...
      Seq.TimeToNextSequence, Seq.tOffset(1) + Seq.CLTime(1));
    Seq.TimeToNextSequence = Seq.tOffset(1) + Seq.CLTime(1);
    Seq.averageRepetitionTime = sum(Seq.tRep(1:end-1)) + Seq.tRelax;
  end

else
  AQ = Seq.AQ;
  TX = Seq.TX;
  Grad = Seq.Grad;
end


if Seq.StartSequence || Seq.PollPPGfast || Seq.GetRawData
  %% Wait for lift movement to finish and store lift position at beginning of sequence
  if Seq.Lift.useLift
    % FIXME: This assumes the previous movement is the same as the currently
    % planned one.
    movementTimeout = HW.Lift.GetMovementTime(Seq.Lift.useSet) + 20e-3;
    err = HW.Lift.WaitForMovement(movementTimeout);
    if err > 0
      error('PD:SampleLift:MovementTime', ...
        'Error %d (%s) waiting for the sample lift movement to finish.', err, err);
    end
  end

  if (isa(HW, 'PD.HWClass') || ~isemptyfield(HW, 'Lift')) && isa(HW.Lift, 'PD.SampleLift') && isvalid(HW.Lift)
    currentLiftPosition = HW.Lift.GetUserPosition(Seq.Lift.checkEncoder);
    currentLiftPositionAbsolute = HW.Lift.GetPosition(Seq.Lift.checkEncoder);
  else
    currentLiftPosition = [];
    currentLiftPositionAbsolute = [];
  end

  if Seq.Lift.useLift
    liftSet = HW.Lift.GetSetInterpreted(Seq.Lift.useSet);
    if ~Seq.Lift.skipLift
      endLiftPosition = currentLiftPosition + HW.Lift.userZero + (2*liftSet.direction-1)*liftSet.distance;
      if (endLiftPosition > HW.Lift.maxPosition) || (endLiftPosition < HW.Lift.minPosition)
        error('PD:SampleLift:ExceedingRange', ...
          'Selected movement would exceed range of sample lift.');
      end
    end
  end

  %% run measurement
  if ~Seq.RawData
    if Seq.Plot || Seq.PlotTR
      [~, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);
    else
      [~, SeqOut, data] = set_sequence(HW, Seq, AQ, TX, Grad);
    end

    iDevice = 1;  % FIXME: Support multiple MMRT devices
    iAQ = find([SeqOut.AQ(:).Device] == iDevice, 1, 'first');  % FIXME: Support multi-channel?

    data(iAQ).lift.position = currentLiftPosition;
    data(iAQ).lift.positionAbsolute = currentLiftPositionAbsolute;
    if SeqOut.Lift.useLift
      data(iAQ).lift.set = liftSet;
    else
      % set field so that lift structures stay compatible
      data(iAQ).lift.set = [];
    end

    if SeqOut.SubtractEmptyMeasurement
      disp('Please remove sample tube and press enter.')
      beep
      if Seq.Lift.useLift && ~Seq.Lift.skipLift
        if isemptyfield(Seq, 'DigitalIO') || isemptyfield(Seq.DigitalIO, 'SetTime')
          Seq.DigitalIO = DigitalIOLift;
        else
          Seq.DigitalIO = add_DigitalIO(Seq.DigitalIO, DigitalIOLift);
        end
      end
      pause
      [~, emptyReference, dataEmpty] = set_sequence(HW, Seq, AQ, TX, Grad);
      SeqOut.emptyReference=emptyReference;
      clear emptyReference;
      emptyReference=SeqOut.emptyReference;
      emptyReference.hParent_Factors=[];
      emptyReference.hParent_fit_exp=[];
      emptyReference.T.hFigure=[];
      emptyReference.T.hParent=[];
      emptyReference.fitExpSettings.hParent=[];
      emptyReference.T.fitExpSettings.hParent=[];
      [pathstr,name,ext] = fileparts(HW.RecoveryCPMG.emptyReferencePath);
      save(fullfile(pathstr, [name, HW.TX(iDevice).CoilName, ext]), 'emptyReference', 'dataEmpty')
      clear emptyReference;
%% fixme
      data.dataEmpty = dataEmpty.data;
      data(iAQ).dataEmpty = dataEmpty.data;
      data(iAQ).data = data(iAQ).data - dataEmpty.data;
    end

  else
    [RawData, SeqOut] = set_sequence(HW, Seq, AQ, TX, Grad);

    if ~iscell(RawData), RawData = {RawData}; end

    iDevice = 1;  % FIXME: Support multiple MMRT devices
    iAQ = find([SeqOut.AQ(:).Device] == iDevice, 1, 'first');  % FIXME: Support multi-channel?

    if SeqOut.SubtractEmptyMeasurement
      disp('Please remove sample tube and press enter.')
      beep
      pause
      [data(iAQ).dataEmpty, ~] = set_sequence(HW, Seq, AQ, TX, Grad);
      RawData{iAQ} = RawData{iAQ} - data(iAQ).dataEmpty;
    end

  end


else
  SeqOut = Seq;
  if isfield(SeqOut, 'data')
    data = SeqOut.data; % FIXME: Where should we really get the data from?
  else
    data = [];
  end

  iDevice = 1;  % FIXME: Support multiple MMRT devices
  iAQ = find([SeqOut.AQ(:).Device] == iDevice, 1, 'first');  % FIXME: Support multi-channel?
end

if SeqOut.PostProcessSequenceLocal
 %% Prepare Data

  if ~SeqOut.RawData
    if (SeqOut.Plot || SeqOut.PlotTR) && ~exist('data_1D', 'var')
      [Seq, data_1D] = get_data_1D(SeqOut, data);
    end
    % Plot
    if SeqOut.Plot
      plot_data_1D(HW, data_1D);
    end

    if SeqOut.PlotTR
      plot_data_1D_TR(HW, data_1D);
    end

    data(iAQ).time_all = data(iAQ).time_all - SeqOut.tAQEchoDelay;
    data(iAQ).time_of_tRep = data(iAQ).time_of_tRep - SeqOut.tAQEchoDelay;

    if ~isempty(SeqOut.SeqAverage) && (SeqOut.SeqAverage.average > 1)
      % nSamples x nAQs x tRep (in Echo train) x nCPMG (+pre-shots)
      newSize = {size(data(iAQ).data,1), size(data(iAQ).data,2), SeqOut.nEcho+2+SeqOut.T2Prepare.tRepsPrepare, []};
      data_CPMG = reshape(data(iAQ).data, newSize{:});
      time_all_CPMG = reshape(data(iAQ).time_all, newSize{:});
      time_of_tRep_CPMG = reshape(data(iAQ).time_of_tRep, newSize{:});

      % nSamples x nAQs x tRep (in Echo train) x averages x nTau1
      newSize = {size(data(iAQ).data,1), size(data(iAQ).data,2), SeqOut.nEcho+2+SeqOut.T2Prepare.tRepsPrepare, SeqOut.SeqAverage.average, []};
      data(iAQ).SeqAverage.data = reshape(data_CPMG(:,:,:,SeqOut.nTau1SteadyState+1:end), newSize{:});
      data(iAQ).SeqAverage.time_all = reshape(time_all_CPMG(:,:,:,SeqOut.nTau1SteadyState+1:end), newSize{:});
      data(iAQ).SeqAverage.time_of_tRep = reshape(time_of_tRep_CPMG(:,:,:,SeqOut.nTau1SteadyState+1:end), newSize{:});

      % mean value of averages (and correction due to phase cycling)
      data(iAQ).data = permute(mean(data(iAQ).SeqAverage.data, 4), [1,2,3,5,4]);
      data(iAQ).dataPhaseCyclingCorrection = 0.5 * ...
        permute(mean(data(iAQ).SeqAverage.data(:,:,:,1:2:end,:), 4) - ...
         mean(data(iAQ).SeqAverage.data(:,:,:,2:2:end,:), 4), [1,2,3,5,4]);
      data(iAQ).time_all = permute(data(iAQ).SeqAverage.time_all(:,:,:,1,:), [1,2,3,5,4]);
      data(iAQ).time_of_tRep = permute(data(iAQ).SeqAverage.time_of_tRep(:,:,:,1,:), [1,2,3,5,4]);
      if ~SeqOut.SeqAverage.SaveSeqAverageData
        data(iAQ).SeqAverage = rmfield(data(iAQ).SeqAverage, {'data','time_all','time_of_tRep'});
      end
    end
     % collect data of all CPMG trains (samples x CPMG x Tau1)
    data(iAQ).SampleEchoTau1Time = reshape(data(iAQ).time_all(:,~isnan(data(iAQ).data(1,:))), size(data(iAQ).data, 1), SeqOut.nEcho - SeqOut.SteadyState_PreShots180 + (SeqOut.tAQFID>=0), []);
    data(iAQ).tStartOftRep = [0 cumsum(SeqOut.tRep(1:end-1).*SeqOut.HW.RX(iDevice).fSample)]./SeqOut.HW.RX(iDevice).fSample;

  else

    data(iAQ).tRep = reshape(find(~isnan(SeqOut.AQ(iAQ).Start(1,:))), [1, max(1, SeqOut.nEcho - SeqOut.SteadyState_PreShots180), SeqOut.SeqAverage.average, max(1, SeqOut.nTau1)]);  % Sample x Echo x Average x Tau1
    data(iAQ).iAQ = ones(size(data(iAQ).tRep));
    if Seq.tAQFID >= 0
      data(iAQ).iAQ(1:(SeqOut.nEcho+1):end) = 2;
    end
    data(iAQ).tStartOftRep = cumsum([0,SeqOut.tRep(1:end-1).*SeqOut.HW.RX(iDevice).fSample],2)./SeqOut.HW.RX(iDevice).fSample;
    data(iAQ).RawTimeOftRep = ...
      ((SeqOut.AQ(iAQ).Start(data(iAQ).iAQ(1),data(iAQ).tRep(1)).*SeqOut.HW.RX(iDevice).fSample.*SeqOut.AQ(iAQ).fSample(data(iAQ).iAQ(1),data(iAQ).tRep(1))) ./ SeqOut.HW.RX(iDevice).fSample ...
       + 0.5 + (0:(SeqOut.AQ(iAQ).nSamples(data(iAQ).iAQ(1),data.tRep(1))-1))).' ...
      ./ SeqOut.AQ(iAQ).fSample(data(iAQ).iAQ(1),data.tRep(1));
    iEchoes = (1+(Seq.tAQFID >= 0)):(size(RawData{iAQ}, 2)/SeqOut.SeqAverage.average/max(1, SeqOut.nTau1));
    data(iAQ).SampleEchoTau1Time(1:size(data(iAQ).RawTimeOftRep,1), iEchoes, 1:max(1, SeqOut.nTau1)) = ...
      reshape(bsxfun(@plus, data(iAQ).RawTimeOftRep.*SeqOut.HW.RX(iDevice).fSample, data(iAQ).tStartOftRep(data(iAQ).tRep(1,:,1,:)) .* SeqOut.HW.RX(iDevice).fSample) ./ SeqOut.HW.RX(iDevice).fSample, ...
              [size(data(iAQ).RawTimeOftRep,1), SeqOut.nEcho, max(1, SeqOut.nTau1)]);
    if Seq.tAQFID >= 0
      data(iAQ).SampleEchoTau1Time(1:SeqOut.AQ(iAQ).nSamples(1,data(iAQ).tRep(1)), 1, 1:max(1, SeqOut.nTau1)) = ...
        bsxfun(@plus, data(iAQ).tStartOftRep(data(iAQ).tRep(1,1,1,:)), ...
        ((SeqOut.AQ(iAQ).Start(1,data(iAQ).tRep(1)).*SeqOut.HW.RX(iDevice).fSample.*SeqOut.AQ(iAQ).fSample(1,data(iAQ).tRep(1))) ./ SeqOut.HW.RX(iDevice).fSample ...
         + 0.5 + (0:(SeqOut.AQ(iAQ).nSamples(1,data(iAQ).tRep(1))-1))).' ...
        ./ SeqOut.AQ(iAQ).fSample(1,data(iAQ).tRep(1)));
      data(iAQ).SampleEchoTau1Time(SeqOut.AQ(iAQ).nSamples(1,data(iAQ).tRep(1))+1:end,1,:) = NaN;
      data(iAQ).SampleEchoTau1Time(SeqOut.AQ(iAQ).nSamples(1,data(iAQ).tRep(2))+1:end,2:end,:) = NaN;
    end

    if ~isempty(SeqOut.SeqAverage) && (SeqOut.SeqAverage.average > 1)
      % phase cycling and averaging result
      data(iAQ).SeqAverage.data = reshape(RawData{iAQ}, ...
        [size(RawData{iAQ},1), 1, SeqOut.nEcho + (SeqOut.tAQFID>=0), SeqOut.SeqAverage.average, max(1, SeqOut.nTau1)]);
      clear RawData
      data(iAQ).data = mean(data(iAQ).SeqAverage.data, 4);
      data(iAQ).dataPhaseCyclingCorrection = 0.5 * ...
        (mean(data(iAQ).SeqAverage.data(:,:,:,1:2:end,:), 4) - ...
         mean(data(iAQ).SeqAverage.data(:,:,:,2:2:end,:), 4));
      if ~SeqOut.SeqAverage.SaveSeqAverageData
        data(iAQ).SeqAverage = rmfield(data(iAQ).SeqAverage, 'data');
      end
    else
      data(iAQ).data = reshape(RawData{iAQ}, ...
        [size(RawData{iAQ},1), 1, SeqOut.nEcho + (SeqOut.tAQFID>=0), 1, max(1, SeqOut.nTau1)]);
      clear RawData
    end
    if Seq.tAQFID >= 0
      data(iAQ).data(SeqOut.AQ.nSamples(1,data(iAQ).tRep(1))+1:end,1,1) = NaN;
      data(iAQ).data(SeqOut.AQ.nSamples(1,data(iAQ).tRep(2))+1:end,1,2:end) = NaN;
      if ~isempty(SeqOut.SeqAverage) && (SeqOut.SeqAverage.average > 1)
        data(iAQ).dataPhaseCyclingCorrection(SeqOut.AQ(iAQ).nSamples(1,data(iAQ).tRep(1))+1:end,1,1) = NaN;
        data(iAQ).dataPhaseCyclingCorrection(SeqOut.AQ(iAQ).nSamples(1,data(iAQ).tRep(2))+1:end,1,2:end) = NaN;
      end
    end

    data(iAQ).lift.position = currentLiftPosition;
    if SeqOut.Lift.useLift
      data(iAQ).lift.set = liftSet;
    else
      % set field so that lift structures stay compatible
      data(iAQ).lift.set = [];
    end
  end
  data(iAQ).SampleEchoTau1 = reshape(data(iAQ).data(:,~isnan(data(iAQ).data(1,:))), ...
    size(data(iAQ).data, 1), SeqOut.nEcho - SeqOut.SteadyState_PreShots180 + (SeqOut.tAQFID>=0), []);
  data(iAQ).tRepsPerCPMG = (SeqOut.T2Prepare.tRepsPrepare+1+SeqOut.nEcho+1);
  data(iAQ).timeOfExcitation = data(iAQ).tStartOftRep(SeqOut.T2Prepare.tRepsPrepare + 1 + ...
    data(iAQ).tRepsPerCPMG*(SeqOut.nTau1SteadyState+(SeqOut.SeqAverage.average*((1:max(1,SeqOut.nTau1))-1))));

  %% post-processing
  iDevice = 1;  % FIXME: Support multiple MMRT devices
  iAQ = find([SeqOut.AQ(:).Device] == iDevice, 1, 'first');  % FIXME: Support multi-channel?

  if ~(SeqOut.StartSequence || SeqOut.PollPPGfast || SeqOut.GetRawData) % only PostProcessSequence, get new raw data
    data(iAQ).SampleEchoTau1 = reshape(data(iAQ).data(:,~isnan(data(iAQ).data(1,:))), ...
      size(data(iAQ).data, 1), SeqOut.nEcho + (SeqOut.tAQFID>=0), []);
  end

  do = true;
  while do
    % use fit_exp to correct the data
    if SeqOut.amplitude2fullScaleReference ~= 1
      data(iAQ).SampleEchoTau1 = data(iAQ).SampleEchoTau1 .* SeqOut.amplitude2fullScaleReference;
    end
    SeqOut.fitExpCorr.Factors = SeqOut.fullScaleReference.Factors(1:min([numel(SeqOut.fullScaleReference.Factors),SeqOut.fullScaleReference.nFactors,SeqOut.nEcho]));
    SeqOut.fitExpCorr.PhaseOffset = SeqOut.fullScaleReference.PhaseOffset;
    if SeqOut.nEcho > 1
      SeqOut.fitExpCorr.EndOffset = 0;
      for iTau1 = 1:max(1, SeqOut.nTau1)
        % correct each CPMG train individually
        [T2_corr, SeqOut.fitExpCorr] = fit_exp(reshape(data(iAQ).SampleEchoTau1(:,:,iTau1), size(data(iAQ).SampleEchoTau1,1), 1, []), ...
          reshape(data(iAQ).SampleEchoTau1Time(:,:,iTau1)-data(iAQ).timeOfExcitation(iTau1), size(data(iAQ).SampleEchoTau1Time,1), 1, []), ...
          SeqOut.fitExpCorr);
        data(iAQ).MeanEchoTau1PhaseCorrection(1,iTau1) = T2_corr.dataPhaseCorrection;
        % reverse correction
        if iTau1 > 1 && all(isnan(T2_corr.dataPhaseCorrected))
          % All data returned by fit_exp could be none if the postprocessing is
          % done before all tau1 steps have been completed. Set the collected
          % data to NaN in that case (irrespective of the dimensions).
          data(iAQ).MeanEchoTau1PhaseCorrected(:,iTau1) = NaN;
          data(iAQ).MeanEchoTau1PhaseCorrectedTime(:,iTau1) = NaN;
        else
          data(iAQ).MeanEchoTau1PhaseCorrected(:,iTau1) = T2_corr.dataPhaseCorrected * exp(1i*T2_corr.dataPhaseCorrection);
          data(iAQ).MeanEchoTau1PhaseCorrectedTime(:,iTau1) = T2_corr.timeCorrected;
        end
      end

      % apply average phase correction
      % data.MeanEchoTau1PhaseCorrected = data.MeanEchoTau1PhaseCorrected * exp(-1i * mean(unwrap(2*data.MeanEchoTau1PhaseCorrection)/2)); % Fixme *2 /2
      % data.MeanEchoTau1PhaseCorrected = data.MeanEchoTau1PhaseCorrected * exp(-1i * mean(unwrap(data.MeanEchoTau1PhaseCorrection)));
      if strcmp(Seq.Recovery, 'Decay')
        % decay
        % use phase of echo train with shortest preparation time as reference for all
        data(iAQ).MeanEchoTau1PhaseCorrected = data(iAQ).MeanEchoTau1PhaseCorrected * exp(-1i * data(iAQ).MeanEchoTau1PhaseCorrection(1,1)); % apply relaxed phase correction
        if real(mean(data(iAQ).MeanEchoTau1PhaseCorrected(1:max(1,floor(end/4)),1))) < 0
          % assume that there is signal at first Tau1
          data(iAQ).MeanEchoTau1PhaseCorrected = -data(iAQ).MeanEchoTau1PhaseCorrected;
          data(iAQ).MeanEchoTau1PhaseCorrection(1,:) = ...
            data(iAQ).MeanEchoTau1PhaseCorrection(1,:) + pi;
        end
      else
        % saturation or inversion
        % use phase of echo train with longest preparation time as reference for all
        iRefPhase = min(iTau1, find(~isnan(data(iAQ).MeanEchoTau1PhaseCorrection), 1, 'last'));
        data(iAQ).MeanEchoTau1PhaseCorrected = data(iAQ).MeanEchoTau1PhaseCorrected * exp(-1i * data(iAQ).MeanEchoTau1PhaseCorrection(1,iRefPhase)); % apply relaxed phase correction
        if real(mean(data(iAQ).MeanEchoTau1PhaseCorrected(1:max(1,floor(end/4)),end))) < 0
          % assume that last Tau1 is relaxed
          data(iAQ).MeanEchoTau1PhaseCorrected = -data(iAQ).MeanEchoTau1PhaseCorrected;
          data(iAQ).MeanEchoTau1PhaseCorrection(1,:) = ...
            data(iAQ).MeanEchoTau1PhaseCorrection(1,:) + pi;
        end
      end
    elseif SeqOut.nTau1 > 1
      % correct the Tau1 "train"
      if isscalar(SeqOut.fitExpCorr.CorrectFrequencyDrift) && SeqOut.fitExpCorr.CorrectFrequencyDrift
        SeqOut.fitExpCorr.CorrectFrequencyDrift = reshape(data(iAQ).SampleEchoTau1Time(:,1,:), size(data(iAQ).SampleEchoTau1Time,1), 1, []);
      end
      % For the phase correction, the time of the excitation pulse is the
      % "reference" time.
      % FIXME: Would it be useful if (some of) the following settings could be
      % overridden by the user?
      SeqOut.fitExpCorr.hasFID = 0;  % If it is *only* an FID, tread it like an echo.
      SeqOut.fitExpCorr.RingFilter = 0;  % no "ring filter" for T1 fit
      SeqOut.fitExpCorr.CorrectFrequencyReferenceTime = 0;  % reference time for FID phase correction is center of excitation pulse
      SeqOut.fitExpCorr.CorrectMeanFrequencyOffset = false;  % T1 measurements are usually "slow". Correct frequency of each AQ window independently.
      [T1_corr, SeqOut.fitExpCorr] = fit_exp(data(iAQ).SampleEchoTau1, ...
        bsxfun(@minus, data(iAQ).SampleEchoTau1Time, reshape(data(iAQ).timeOfExcitation, 1, 1, [])), ...
        SeqOut.fitExpCorr);

      data(iAQ).MeanEchoTau1PhaseCorrected(1,:) = T1_corr.dataPhaseCorrected;
      if real(data(iAQ).MeanEchoTau1PhaseCorrected(1,end)) < 0
        % assume that last Tau1 is relaxed
        data(iAQ).MeanEchoTau1PhaseCorrected = -data(iAQ).MeanEchoTau1PhaseCorrected;
      end
      data(iAQ).MeanEchoTau1PhaseCorrectedTime(1,:) = T1_corr.timeCorrected;
    else
      % only one single Echo (don't do any corrections)
      data(iAQ).MeanEchoTau1PhaseCorrectedTime = SeqOut.tEcho;
      data(iAQ).MeanEchoTau1PhaseCorrected = mean(data(iAQ).SampleEchoTau1);
    end


    % data.MeanEchoTau1 = reshape(mean(data.SampleEchoTau1, 1), size(data.SampleEchoTau1,2), []);
    % % data.MeanEchoTau1PhaseOffset = mean(angle(data.MeanEchoTau1(1:min(5,end),1)));              % Phase at first Tau1
    % data.MeanEchoTau1PhaseOffset = mean(angle(data.MeanEchoTau1(1:ceil(min(4,end/8)),end)));  % Phase at long Tau1
    % data.MeanEchoTau1PhaseCorrected_old = data.MeanEchoTau1 * exp(-1i*data.MeanEchoTau1PhaseOffset);
    % data.EchoTime = cumsum(repmat(SeqOut.tEcho, [SeqOut.nEcho, max(1,SeqOut.nTau1)]), 1);
    % data.Tau1Time = repmat(SeqOut.Tau1(:).', SeqOut.nEcho, 1);

    data(iAQ).EchoTime = data(iAQ).MeanEchoTau1PhaseCorrectedTime;
    if SeqOut.nTau1 > 0
      % For the T1 fit, the time between inversion (preparation) pulse and
      % excitation pulse is the "reference" time.
      data(iAQ).Tau1Time = repmat(SeqOut.Tau1(:).', size(data(iAQ).EchoTime, 1), 1);
    end
    do = SeqOut.MeasureFullScaleReference && SeqOut.amplitude2fullScaleReference == 1;
    %% get correction factors for n-th Echo
    if do
      FSR = SeqOut.fullScaleReference;
      FSR.RefEchos = data(iAQ).MeanEchoTau1PhaseCorrected(:,FSR.SelectedTau1);
      FSR.time = data(iAQ).EchoTime(:,FSR.SelectedTau1(1));
      FSR.RefEchosMean = mean(FSR.RefEchos, 2);
      FSR.RefEchosFitEndIndex = min([find(FSR.RefEchosMean<0, 1, 'first'), numel(FSR.RefEchosMean)]);
      % [FSR.T, FSR.fitExpSettings]=fit_exp(abs(FSR.RefEchosMean(FSR.fitStart:FSR.RefEchosFitEndIndex)),FSR.time(FSR.fitStart:FSR.RefEchosFitEndIndex),FSR.hParent_fit_exp,0,0,0,0);
      [FSR.T, FSR.fitExpSettings] = fit_exp(FSR.RefEchosMean(FSR.fitStart:FSR.RefEchosFitEndIndex), ...
        FSR.time(FSR.fitStart:FSR.RefEchosFitEndIndex), FSR.hParent_fit_exp, 0, 1, 0, 1);
      FSR.bestFitDecaySingle = bsxfun(@plus, ...
        FSR.T.xminSingle(1), ...
        bsxfun(@times, FSR.T.xminSingle(2), exp(-bsxfun(@rdivide, FSR.time, FSR.T.xminSingle(3)))));
      FSR.PhaseOffset = FSR.T.dataPhaseCorrection;
      FSR.Factors = FSR.bestFitDecaySingle ./ (FSR.RefEchosMean .* exp(-1i*FSR.PhaseOffset));  % rel amp
      FSR.FactorsStd = std(bsxfun(@rdivide, FSR.bestFitDecaySingle, FSR.RefEchos), 0, 2);  % rel amp
      % scale amplitude to reference (100% water)
      FSR.Amplitude = FSR.T.xminSingle(2);          % reference amplitude for water sample

      if FSR.hParent_Factors > 0
        hf = figure(FSR.hParent_Factors);
        clf(hf);
        hax(1) = subplot(3,1,1, 'Parent', hf);
        plot(hax(1), abs(FSR.Factors));
        ylabel(hax(1), {'Correction Factor to' , 'Single Exponential Fit'});
        ylim(hax(1), [0.8, 1.2]);
        grid(hax(1), 'on');
        title(hax(1), ['FSR.Amplitude = ' num2str(FSR.Amplitude*1e9,6) ' nT']);

        hax(2)=subplot(3,1,2,'parent',hf);
        plot(hax(2), FSR.FactorsStd);
        ylabel(hax(2), {'Correction Factor STD to' , 'Single Exponential Fit'});
        grid(hax(2), 'on');

        hax(3)=subplot(3,1,3, 'Parent', hf);
        plot(hax(3), angle(FSR.Factors));
        xlabel(hax(3), 'Number of Echo');
        ylabel(hax(3), {'Correction Phase to' , 'Single Exponential Fit'});
        title(hax(3), ['FSR.PhaseOffset = ' num2str(FSR.PhaseOffset) ' rad'])
        grid(hax(3), 'on');
        linkaxes(hax(1:3), 'x')
      end
      SeqOut.fullScaleReference = FSR;
      clear FSR;
      fullScaleReference = SeqOut.fullScaleReference;
      fullScaleReference.hParent_Factors = [];
      fullScaleReference.hParent_fit_exp = [];
      fullScaleReference.T.hFigure = [];
      fullScaleReference.T.hParent = [];
      fullScaleReference.fitExpSettings.hParent = [];
      fullScaleReference.T.fitExpSettings.hParent = [];
      [pathstr,name,ext] = fileparts(HW.RecoveryCPMG.fullScaleReferencePath);
      save(fullfile(pathstr, [name, HW.TX(iDevice).CoilName, ext]), 'fullScaleReference');
      clear fullScaleReference;
      SeqOut.amplitude2fullScaleReference = 1 / SeqOut.fullScaleReference.Amplitude ...
        * SeqOut.fullScaleReference.sampleCrossSection * min(SeqOut.fullScaleReference.thicknessSlice, abs(diff(HW.Grad(iDevice).ImageVol(3:4)))) ...
        / SeqOut.sampleCrossSection / min(SeqOut.thicknessSlice, abs(diff(HW.Grad(iDevice).ImageVol(3:4))));
    end
  end

  if SeqOut.PlotT1T2 > 0 && SeqOut.nTau1 > 1 && SeqOut.nEcho > 1
    %% plot corrected measurement amps as pseudocolor plot (2d)
    hf80 = figure(SeqOut.PlotT1T2);
    hf80 = clf(hf80);
    axf80 = axes('Parent', hf80);
    % imagesc(data.EchoTime(:,1),data.Tau1Time(1,:),real(data.MeanEchoTau1PhaseCorrected))
    hp = pcolor(axf80, data(iAQ).EchoTime, data(iAQ).Tau1Time, real(data(iAQ).MeanEchoTau1PhaseCorrected));
    set(hp, 'ZData', real(data(iAQ).MeanEchoTau1PhaseCorrected));
    set(axf80, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'linear');
    ylabel(axf80, 'Tau1 / sec');
    xlabel(axf80, 'Echo time / sec');
    shading(axf80, 'interp');
    colorbar('peer', axf80);

    % figure(6)
    %imagesc(data.EchoTime(:,1),data.Tau1Time(1,:),real(data.MeanEchoTau1PhaseCorrected))
    % % pcolor(data.EchoTime,data.Tau1Time,real(data.MeanEchoTau1PhaseCorrected))
    % xlabel('Tau1 / sec')
    % ylabel('Echo time / sec')
    % shading interp
    % colorbar
    %
    % pause(0.5)

    % plot corrected measurement amps as surface (3d)
    hf81 = figure(SeqOut.PlotT1T2+1);
    hf81 = clf(hf81);
    axf81 = axes('Parent', hf81);
    surface(data(iAQ).EchoTime, data(iAQ).Tau1Time, real(data(iAQ).MeanEchoTau1PhaseCorrected), 'Parent', axf81);
    set(axf81, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'linear');
    ylabel(axf81, 'Tau1 / sec');
    xlabel(axf81, 'Echo time / sec');
    zlabel(axf81, 'real part of Echo amplitude')
    view(axf81, 3); % 3d
    colorbar('peer', axf81);
    grid(axf81, 'on');
    shading(axf81, 'interp');

    % plot residual amps as surface (3d)
    hf81 = figure(SeqOut.PlotT1T2+101);
    hf81 = clf(hf81);
    axf81 = axes('Parent', hf81);
    surface(data(iAQ).EchoTime, data(iAQ).Tau1Time, imag(data(iAQ).MeanEchoTau1PhaseCorrected), 'Parent', axf81);
    set(axf81, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'linear');
    ylabel(axf81, 'Tau1 / sec');
    xlabel(axf81, 'Echo time / sec');
    zlabel(axf81, 'imag part of Echo amplitude (residuals)')
    view(axf81, 3); % 3d
    colorbar('peer', axf81);
    grid(axf81, 'on');
    shading(axf81, 'interp');
  end

  %% T1 fit to first (or last) echoes
  if SeqOut.FitT1
    if SeqOut.nTau1 > 2
      % [T1, SeqOut.fitExpT1] = fit_exp(data.MeanEchoTau1PhaseCorrected(1,:).*data.MeanEchoTau1PhaseCorrection, data.Tau1Time(1,:), fitExpSettings);
      [T1, SeqOut.fitExpT1] = fit_exp(data(iAQ).MeanEchoTau1PhaseCorrected(1,:), data(iAQ).Tau1Time(1,:), SeqOut.fitExpT1);
      if ishghandle(T1.hParent, 'figure')
        set(T1.hParent, 'Name', 'T1 Fit');
      end
      % t1 = T1.tau;
      T1.half = log(2)*T1.tau;
      T1.half1 = log(2)*T1.tau1;
      T1.half2 = log(2)*T1.tau2;
      data(iAQ).T1 = T1;
    else
      disp('nTau1 too low')
    end
  end

  %% T2 to selected CPMG trains
  if SeqOut.FitT2 && numel(SeqOut.FitT2AtTau1) > 0
    if SeqOut.nEcho > 3
      if SeqOut.PlotT2 > 0
        fh = figure(SeqOut.PlotT2);
        clf(fh)
        set(fh, 'Name', 'T2 Fit');
        ax(1) = axes('Parent', fh);
        hold(ax(1), 'on');
      end
      if isemptyfield(SeqOut.fitExpT2, 'omitFirstnEchoes')
        SeqOut.fitExpT2.omitFirstnEchoes = SeqOut.fitExpCorr.omitFirstnEchoes;
      end
      for t = SeqOut.FitT2AtTau1
        % fitting exponential functions to the data
        if SeqOut.PlotT2 > 0
          SeqOut.fitExpT2.hParent = figure(SeqOut.PlotT2*10+t);
        else
          SeqOut.fitExpT2.hParent = 0;
        end
        [T2(t), SeqOut.fitExpT2] = fit_exp(data(iAQ).MeanEchoTau1PhaseCorrected(:,t), data(iAQ).EchoTime(:,t), SeqOut.fitExpT2);
        if SeqOut.PlotT2 > 0
          set(SeqOut.fitExpT2.hParent, 'Name', sprintf('T2 Fit at tau_1 #%d', t));
          plot(ax(1), T2(t).timeCorrected, T2(t).dataPhaseCorrectedReal, 'b');
          legendEntry = {'real'};
          if SeqOut.fitExpT2.DoubleExp
            plot(ax(1), T2(t).functionTime, T2(t).functionAmpDouble, 'b-.');
            legendEntry{2} = 'fit double';
          elseif SeqOut.fitExpT2.SingleExp
            plot(ax(1), T2(t).functionTime, T2(t).functionAmpSingle, 'b--');
            legendEntry{2} = 'fit single';
          end
        end
      end
      if SeqOut.PlotT2 > 0
        if numel(SeqOut.FitT2AtTau1) > 0
          % legend(ax(1),'real');
          xlabel(ax(1), 'time / s');
          ylabel(ax(1), 'amplitude');
          title(ax(1), {[' ', ' ', ' ']});
          grid(ax(1), 'on');
          legend(ax(1), legendEntry);
        end
        hold(ax(1), 'off');
      end

      data(iAQ).T2 = T2;
    else
      disp('nEcho too low')
    end
  end

  % SeqOut.iLaplace2D.T1Start=SeqOut.Tau1Start/1;
  % SeqOut.iLaplace2D.T1End=SeqOut.Tau1End*5;
  % SeqOut.iLaplace2D.T2Start=SeqOut.tEcho*1;
  % SeqOut.iLaplace2D.T2End=SeqOut.tEcho*SeqOut.nEcho*5;
  % SeqOut.iLaplace2D.n_Lcomp=min(SeqOut.nTau1,30);


end

end


function Seq = SetPowerSupplyAmp(Seq, nn, oldPreStartPPGFcn, HW, gradAmp)
%% Function that configures the power supply for the calculated amplitude

if ~isempty(oldPreStartPPGFcn)
  % call previous PreStartPPGFcn
  Seq = oldPreStartPPGFcn(Seq, nn);
end

if nn ~= 1
  % only send commands on first average loop
  return;
end

iDevice = 1;  % FIXME: Support pulse program at secondary devices.

oldLocked = HW.PowerSupply.isLocked();
if ~oldLocked
  % Set lock or changing the limit for the current doesn't have an effect.
  HW.PowerSupply.SetLock(true);
  % FIXME: Do we need to reset to the previous lock state after the measurement?
end

% Compensate reduced coil efficiency such that the actual current output is
% limitted correctly digitally
if HW.Grad(iDevice).PaCurrentControlled(HW.Grad(iDevice).Slice.channel)
  HW.PowerSupply.PaIoutMax = gradAmp * HW.Grad(iDevice).Amp2PaIout(HW.Grad(iDevice).Slice.channel) ...
    / HW.PowerSupply.AnalogOverDriveFactor;
  % HW.PowerSupply.PaIoutMax = gradAmp * HW.Grad(iDevice).Amp2LoadIin(HW.Grad(iDevice).Slice.channel) ...
  %   / HW.PowerSupply.AnalogOverDriveFactor;
else
  error('PD:SetPowerSupplyAmp:NotCC', ...
    'Gradient amplifier must be current controlled.');
end

% FIXME: Do we need this for the remote pin to work?
HW.PowerSupply.SetLock(false);

end
