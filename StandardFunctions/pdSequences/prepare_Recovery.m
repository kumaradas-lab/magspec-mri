function [HW, Seq, AQ, TX, Grad] = prepare_Recovery(HW, Seq, AQ, TX, Grad)
%% Prepare sequence using a inversion or saturation method to measure the T1 time
%
%     [HW, Seq, AQ, TX, Grad] = prepare_Recovery(HW, Seq, AQ, TX, Grad)
%
%
% INPUT:
%
%   HW
%       HW object or structure.
%
%   Seq
%       Structure with the following (optional) fields. If these fields are
%       omitted or empty, default values might be used.
%
%     nPreLoops
%         Number of pre-shots without acquisition.
%         (Default: 1)
%
%     nFlips
%         Number of different recovery times after preparation pulse.
%         (Default: Seq.Loops-Seq.nPreLoops)
%
%     Recovery
%         Type of the recovery experiment. Must be either 'Saturation' or
%         'Inversion'.
%         (Default: 'Inversion')
%
%     inversionPulse
%         Pulse shape function of the inversion pulse (if used).
%         (Default: @Pulse_Rect_Composite180)
%
%     saturationPulse
%         Pulse shape function of the saturation pulse (if used).
%         (Default: @Pulse_Rect_Composite90)
%
%     preparationPulseProperties
%         Additional properties that should be passed to the pulse shape
%         function.
%
%     nSaturationPulses
%         Number of saturation pulses (if used).
%         (Default: 1)
%
%     tFlip
%         Vector with durations of recovery time in seconds (i.e., the time
%         between the center of the preparation pulse and the center of the
%         excitation pulse). If it is a scalar, a vector with "reasonable" times
%         is automatically created (see below).
%         (Default: 10e-3)
%
%     tFlipLog
%         Boolean value to indicate if the vector with recovery times should be
%         created with logarithmically increasing steps (true) or equally spaced
%         steps (false).
%         (Default: false)
%
%     tFlipStart
%         If tFlip is a scalar, this is the shortest recovery time in seconds
%         for the vector of recovery times.
%         (Default: Seq.tFlip)
%
%     tFlipEnd
%         If tFlip is a scalar, this is the longest recovery time in seconds
%         for the vector of recovery times.
%         (Default: Seq.tFlip*Seq.nFlips)
%
%     tShiftPreparation
%         Optionally, shift all preparation pulses by a fixed duration in
%         seconds. This value adds to the recovery time.
%         (Default: 0)
%
%     tRelax
%         Relaxation time in seconds between experiments with increments of the
%         recovery time.
%         (Default: Seq.LoopsBreak)
%
%     tSaturationGrad
%         Duration in seconds of a spoiler pulse in the saturation recovery
%         experiment.
%         (Default: min(20e-3,
%         Seq.tFlip(1)-Seq.tGradEC-HW.Grad(iDevice).tRamp*2-HW.Grad(iDevice).tEC*2-HW.tFlip90Def-0.5e-3)
%         )
%
%     AmpSaturationGrad
%         Amplitude in Tesla of the spoiler pulse in the saturation recovery
%         experiment. A gradient pulse with that amplitude is emitted in x-, y-,
%         *and* z-direction.
%         (Default: min(100e-3, min(HW.Grad(iDevice).MaxAmp(1:3))/5) )
%
%     useLateSpoiler
%         Boolean value that selects whether the spoiler pulse in the saturation
%         recovery experiment should be emitted Seq.tGradEC seconds after the
%         preparation pulse (false), or immediately before the excitation pulse
%         (true).
%         (Default: false)
%
%     InversionSpoilSize
%         Size in meters of the spoiler after the preparation pulse in the
%         inversion recovery experiment.
%         (Default: 1e-3)
%
%     InversionSpoilLength
%         Duration in seconds of the spoiler after the preparation pulse in the
%         inversion recovery experiment.
%         (Default: 2e-3)
%
%   AQ
%       Structure with acquisition properties.
%
%   TX
%       Structure with rf pulse properties.
%
%   Grad
%       Structure with gradient properties.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% if isfield(Seq, 'LoopName') && isfield(Seq, 'LoopNameCount') && ...
%     ~strcmp(Seq.LoopName{Seq.LoopNameCount}, 'normal')
%   % only change normal loops when called from sequence_Spin_Echo
%   return;
% end
if ~strcmp(Seq.LoopName{Seq.LoopNameCount}, 'normal')
  % only change normal loops
  return;
end

%% Default parameters
Seq = set_EmptyField(Seq, 'nPreLoops', 1);
Seq = set_EmptyField(Seq, 'nFlips', Seq.Loops-Seq.nPreLoops);
Seq = set_EmptyField(Seq, 'inversionPulse', @Pulse_Rect_Composite180);
Seq = set_EmptyField(Seq, 'saturationPulse', @Pulse_Rect_Composite90);
if isemptyfield(Seq, 'preparationPulseProperties')
  Seq.preparationPulseProperties = struct();
end
Seq = set_EmptyField(Seq, 'nSaturationPulses', 1);
Seq = set_EmptyField(Seq, 'tFlip', 10e-3);
Seq = set_EmptyField(Seq, 'tFlipLog', 0);
if isscalar(Seq.tFlip)
  Seq = set_EmptyField(Seq, 'tFlipStart', Seq.tFlip);
  Seq = set_EmptyField(Seq, 'tFlipEnd', Seq.tFlip*Seq.nFlips);
  % Create vector with times for the inversion pulse
  if Seq.tFlipLog
    Seq.tFlip = logspace(log10(Seq.tFlipStart), log10(Seq.tFlipEnd), Seq.nFlips);
  else
    Seq.tFlip = linspace(Seq.tFlipStart, Seq.tFlipEnd, Seq.nFlips);
  end
end
if isemptyfield(Seq, 'tShiftPreparation')
  Seq.tShiftPreparation = 0;
end
Seq = set_EmptyField(Seq, 'tRelax', Seq.LoopsBreak);  % FIXME: Is this unused?
Seq = set_EmptyField(Seq, 'Recovery', 'Inversion');
Seq = set_EmptyField(Seq, 'tGradEC', 5e-3);
iDevice = Seq.AQSlice(1).iDevice;
Seq = set_EmptyField(Seq, 'tSaturationGrad', ...
  min(20e-3, Seq.tFlip(1)-Seq.tGradEC-HW.Grad(iDevice).tRamp*2-HW.Grad(iDevice).tEC*2-HW.tFlip90Def-0.5e-3));
Seq = set_EmptyField(Seq, 'AmpSaturationGrad', ...
  min(100e-3, min(HW.Grad(iDevice).MaxAmp(1:3))/5));
if isemptyfield(Seq, 'useLateSpoiler')
  Seq.useLateSpoiler = false;
end


Seq.tFlip = Seq.tFlip([ones(1,Seq.nPreLoops),1:end]);  % repeat the first flip time for the pre-loops

% keep timing between end of Turbo blocks and preparation pulse of
% subsequent turbo block constant for each preparation time
Seq.AQSlice(1).TurboBreak = Seq.tRelax + Seq.tFlip(Seq.Loop);
Seq.tRep(Seq.P90tReps(2:end)-1) = Seq.tEcho + Seq.AQSlice(1).TurboBreak;
Seq.LoopsBreak = Seq.AQSlice(1).TurboBreak;

%% Add pulses
switch Seq.Recovery
  case 'Inversion'
    Seq.iLaplace1D.Problem = 'Inversion';

    % default parameters specific to inversion recovery
    if isemptyfield(Seq, 'InversionSpoilSize')
      % Size of spoiler at first inversion pulse in meter
      Seq.InversionSpoilSize = 1e-3;
    end
    if isemptyfield(Seq, 'InversionSpoilLength')
      % Length of spoiler at first inversion pulse in seconds
      Seq.InversionSpoilLength = 2e-3;
    end

    tSlice = numel(Seq.Slice);
    % 90 degrees slice excitation
    % Set the pulse shape for the initial inversion pulse:
    Seq.Slice(tSlice+1).Pulse.Function = Seq.inversionPulse;
    Seq.Slice(tSlice+1).Pulse.FlipAngleComposite = Seq.AQSlice(1).inversionFlipAngleComposite;
    Seq.Slice(tSlice+1).Pulse.FlipPhaseComposite = Seq.AQSlice(1).inversionFlipPhaseComposite;
    % Make it an inversion pulse:
    Seq.Slice(tSlice+1).Pulse.FlipAngle = 180;
    % Basically without slice selection
    Seq.Slice(tSlice+1).Thickness = 1e12;
    % Set the center of the inversion pulse to "tFlip" seconds prior to the the
    % center of the excitation pulse. "Seq.Loop" is the loop counter. Thus, use
    % a higher tFlip every loop.
    Seq.Slice(tSlice+1).CenterOfPulse = Seq.tEcho/2 - Seq.tFlip(Seq.Loop) - Seq.tShiftPreparation;
    % Add this pulse to the tRep with the excitation pulse(s):
    Seq.Slice(tSlice+1).UseAtRepetitionTime = Seq.P90tReps;
    preparationPulseProperties = fieldnames(Seq.preparationPulseProperties);
    for iProp = 1:numel(preparationPulseProperties)
      % copy additional properties for pulse shape function
      Seq.Slice(tSlice+1).Pulse.(preparationPulseProperties{iProp}) = ...
        Seq.preparationPulseProperties.(preparationPulseProperties{iProp});
    end

    % Seq.Slice(tSlice+1).UseCoordinate=Seq.AQSlice(1).SliceCoordinate;
    % Seq.Slice(tSlice+1).GradTimeDelay=Seq.AQSlice(1).SliceTimeDelay;
    % Seq.Slice(tSlice+1).Pulse.Phase=Seq.tRepTurboBlockPhaseExcitation(Seq.Slice(1).UseAtRepetitionTime)+Seq.AQSlice(1).excitationPhase;
    % Seq.Slice(tSlice+1).Pulse.PhaseIncrement=Seq.AQSlice(1).excitationPhaseIncrement;
    % Seq.Slice(tSlice+1).CenterOfRephase=Seq.Read(1).CenterOfRephase-Seq.tReadoutOffset-mean(Seq.Read(1).GradTimeDelayOffset);
    % Seq.Slice(tSlice+1).GradRephaseLength=Seq.Read(1).GradRephaseLength;
    % Seq.Slice(tSlice+1).MaxGradAmp=Seq.MaxGradAmpSlice;
    % Seq.Slice(tSlice+1).distance=Seq.AQSlice(1).Center2OriginImage(1);
    % Seq.Slice(tSlice+1).GradTimeIntegralRephaseOffset=Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset;
    % Seq.Slice(tSlice+1).GradSign=Seq.AQSlice(1).SliceGradSign*Seq.AQSlice(1).GradSign;
    % Seq.Slice(tSlice+1).GradDephaseSign=-Seq.AQSlice(1).SliceGradSign*Seq.AQSlice(1).GradSign;
    % Seq.Slice(tSlice+1).GradRephaseSign=-Seq.AQSlice(1).SliceGradSign*Seq.AQSlice(1).GradSign;

    % generate slice parameters
    Seq = get_SliceParameter(Seq, HW);


    % FIXME: Add support for Seq.useLateSpoiler
    tPhase = numel(Seq.Phase);
    % 180 degrees spoiler encoding aligned with read rephase
    % Size of the spoiler in meter:
    Seq.Phase(tPhase+1).sizePhase = Seq.InversionSpoilSize;
    % The spoiler does not actually encode the image in phase direction:
    Seq.Phase(tPhase+1).nPhase = 1;
    % "Default" settings for a spoiler
    Seq.Phase(tPhase+1).PhaseOS = 2;
    Seq.Phase(tPhase+1).StepOrder = [1,1];
    % Length of the spoiler in seconds:
    Seq.Phase(tPhase+1).GradDephaseLength = Seq.InversionSpoilLength;
    % Apply spoiler in read direction (1):
    Seq.Phase(tPhase+1).UseCoordinate = 1;
    % Add the spoiler in the same tRep as the inversion pulse:
    Seq.Phase(tPhase+1).UseAtRepetitionTime = Seq.Slice(tSlice+1).UseAtRepetitionTime;
    % Sign of the spoiler gradient amplitude (arbitrary):
    Seq.Phase(tPhase+1).GradDephaseSign = 1;
    Seq.Phase(tPhase+1).CenterOfDephase = Seq.Slice(tSlice+1).CenterOfPulse + Seq.Slice(tSlice+1).GradLength/2 + Seq.Phase(tPhase+1).GradDephaseLength/2;
    % A rephase gradient before the inversion pulse is not necessary and it
    % won't be added to the sequence below. But add it for completeness anyway:
    Seq.Phase(tPhase+1).GradRephaseLength = Seq.Phase(tPhase+1).GradDephaseLength;
    Seq.Phase(tPhase+1).GradRephaseSign = Seq.Phase(tPhase+1).GradDephaseSign;
    Seq.Phase(tPhase+1).CenterOfRephase = Seq.Slice(tSlice+1).CenterOfPulse - Seq.Slice(tSlice+1).GradLength/2 - Seq.Phase(tPhase+1).GradRephaseLength/2;

    % generate phase parameters
    Seq = get_PhaseParameter(Seq, HW);

    TX = add_TX(TX,Seq.Slice(tSlice+1).TX);  % Add the inversion pulse
    Grad = add_Grad(Grad, Seq.Slice(tSlice+1).Grad);  % Add the slice selection gradient
    Grad = add_Grad(Grad, Seq.Phase(tPhase+1).GradDephase);  % Add the spoiler

  case 'Saturation'
    Seq.iLaplace1D.Problem = 'Saturation';
    Seq.tSaturat = HW.tFlip90Def * Seq.saturationPulse(HW, 'Amp');

    % FIXME: Use get_SliceParameter to create preparation pulse(s).
    Seq.PulsePreparation.Function = Seq.saturationPulse;
    % Make it an saturation pulse:
    Seq.PulsePreparation.FlipAngle = pi * HW.tFlip90Def ...
      / (HW.TX(iDevice).Amp2FlipPiIn1Sec / HW.TX(iDevice).AmpDef);
    % % Basically without slice selection
    % Seq.Slice(tSlice+1).Thickness = 1e12;
    % Set the center of the inversion pulse to "tFlip" seconds prior to the the
    % center of the excitation pulse. "Seq.Loop" is the loop counter. Thus, use
    % a higher tFlip every loop.
    Seq.PulsePreparation.CenterOfPulse = Seq.tEcho/2 - Seq.tFlip(Seq.Loop);
    % % Add this pulse to the tRep with the excitation pulse(s):
    % Seq.Slice(tSlice+1).UseAtRepetitionTime = Seq.P90tReps;
    Seq.PulsePreparation.Bandwidth = 1/Seq.tSaturat*Seq.saturationPulse(HW,'Time');
    Seq.PulsePreparation.MaxNumberOfSegments = 51;
    Seq.PulsePreparation.MaxLength = Seq.tSaturat;
    % only applies to spin-locking pulses
    if isemptyfield(Seq.PulsePreparation, 'DurationSpinLock')
      Seq.PulsePreparation.DurationSpinLock = Seq.tFlip(Seq.Loop);
    end
    preparationPulseProperties = fieldnames(Seq.preparationPulseProperties);
    for iProp = 1:numel(preparationPulseProperties)
      % copy additional properties for pulse shape function
      Seq.PulsePreparation.(preparationPulseProperties{iProp}) = ...
        Seq.preparationPulseProperties.(preparationPulseProperties{iProp});
    end

    pulsePrepare = Seq.saturationPulse(HW, 0, Seq.PulsePreparation);

    % % generate slice parameters
    % Seq = get_SliceParameter(Seq, HW);

    % Gradients
    if Seq.tSaturat/2+HW.Grad(iDevice).tEC+Seq.tSaturationGrad+HW.Grad(iDevice).tRamp+HW.Grad(iDevice).tEC+Seq.tSaturat/2+50e-6 > Seq.tFlip(1)
      error('Increase Seq.tFlipStart or decrease Seq.tSaturationGrad');
    end
    temp = nan(size(Seq.tRep));
    temp(Seq.P90tReps)=1;
    for t = 1:HW.Grad(iDevice).n
      if t <= 3
        if Seq.useLateSpoiler
          Seq.GradRecovery(t).Time = ...
            cumsum([0; ...
            HW.Grad(iDevice).tRamp; ...
            Seq.tSaturationGrad-HW.Grad(iDevice).tRamp*1; ...
            HW.Grad(iDevice).tRamp]) * temp ...
            + Seq.Slice(1).CenterOfPulse - Seq.Slice(1).GradLength/2 - HW.Grad.tEC ...
            - (Seq.tSaturationGrad + HW.Grad(iDevice).tRamp) - 1e-3;
        else
          Seq.GradRecovery(t).Time = ...
            cumsum([Seq.tSaturat/2+HW.Grad(iDevice).tEC; ...
            HW.Grad(iDevice).tRamp; ...
            Seq.tSaturationGrad-HW.Grad(iDevice).tRamp*1; ...
            HW.Grad(iDevice).tRamp]) * temp ...
            - Seq.tFlip(Seq.Loop) + Seq.Slice(1).CenterOfPulse;
        end
        Seq.GradRecovery(t).Amp = Seq.AmpSaturationGrad * [0;1;1;0;] * temp;
      else
        Seq.GradRecovery(t).Time = NaN;
        Seq.GradRecovery(t).Amp = NaN;
      end
    end

    if Seq.GradRecovery(1).Time(4,1) > Seq.Slice(1).CenterOfPulse-Seq.Slice(1).GradLength/2-HW.Grad.tEC
      error('prepare_Recovery:tGradEC', 'Saturation gradient too close to Pulse (Seq.tGradEC)')
    end

    for tt = 2:Seq.nSaturationPulses
      if tt==2, pulsePrepareTemp = pulsePrepare; end
      pulsePrepare.Start     =cat(1, ...
        pulsePrepare.Start, ...
        pulsePrepareTemp.Start-(tt-1)*(Seq.tSaturationGrad+Seq.tSaturat*2+HW.Grad(iDevice).tEC*2+HW.Grad(iDevice).tRamp+Seq.tGradEC));
      pulsePrepare.Duration  =cat(1, pulsePrepare.Duration, pulsePrepareTemp.Duration);
      pulsePrepare.Amplitude =cat(1, pulsePrepare.Amplitude, pulsePrepareTemp.Amplitude);
      pulsePrepare.Frequency =cat(1, pulsePrepare.Frequency, pulsePrepareTemp.Frequency);
      pulsePrepare.Phase     =cat(1, pulsePrepare.Phase, pulsePrepareTemp.Phase);
      if tt == 2, GradRecoveryTemp = Seq.GradRecovery;end
      for t=1:HW.Grad(iDevice).n
        Seq.GradRecovery(t).Time = cat(1, ...
          Seq.GradRecovery(t).Time, ...
          GradRecoveryTemp(t).Time-(tt-1)*(Seq.tSaturationGrad+Seq.tSaturat*2+HW.Grad(iDevice).tEC*2+HW.Grad(iDevice).tRamp+Seq.tGradEC));
        % Seq.GradRecovery(t).Amp=cat(1,Seq.GradRecovery(t).Amp,GradRecoveryTemp(t).Amp*(((Seq.nSaturationPulses-tt+1))/Seq.nSaturationPulses).^2);
        Seq.GradRecovery(t).Amp = cat(1, Seq.GradRecovery(t).Amp, GradRecoveryTemp(t).Amp);
      end
    end

    Seq.TXRecovery.Start = NaN(size(pulsePrepare.Start,1), size(Seq.tRep,2));
    Seq.TXRecovery.Duration = Seq.TXRecovery.Start;
    Seq.TXRecovery.Amplitude = Seq.TXRecovery.Start;
    Seq.TXRecovery.Frequency = Seq.TXRecovery.Start;
    Seq.TXRecovery.Phase = Seq.TXRecovery.Start;

    Seq.TXRecovery.Duration(1:size(pulsePrepare.Start,1),Seq.P90tReps) = ...
      repmat(pulsePrepare.Duration, 1, numel(Seq.P90tReps));
    Seq.TXRecovery.Start(1:size(pulsePrepare.Start,1),Seq.P90tReps) = ...
      repmat(pulsePrepare.Start-Seq.tFlip(Seq.Loop)+Seq.Slice(1).CenterOfPulse, 1, numel(Seq.P90tReps)) ...
      - Seq.tShiftPreparation;
    Seq.TXRecovery.Amplitude(1:size(pulsePrepare.Start,1),Seq.P90tReps) = ...
      repmat(pulsePrepare.Amplitude, 1, numel(Seq.P90tReps));
    Seq.TXRecovery.Frequency(1:size(pulsePrepare.Start,1),Seq.P90tReps) = ...
      repmat(pulsePrepare.Frequency, 1, numel(Seq.P90tReps));
    Seq.TXRecovery.Phase(1:size(pulsePrepare.Start,1),Seq.P90tReps) = ...
      repmat(pulsePrepare.Phase, 1, numel(Seq.P90tReps));

    TX = add_TX(TX, Seq.TXRecovery);
    Grad = add_Grad(Grad, Seq.GradRecovery);

  otherwise
    error('Seq.Recovery must be either "Inversion" or "Saturation"!');

end

end
