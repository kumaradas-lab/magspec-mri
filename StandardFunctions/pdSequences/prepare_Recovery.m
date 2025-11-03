function [HW, Seq, AQ, TX, Grad] = prepare_Recovery(HW, Seq, AQ, TX, Grad)
%% Prepare sequence using a inversion or saturation method to measure the T1 time
%
%     [HW, Seq, AQ, TX, Grad] = prepare_Recovery(HW, Seq, AQ, TX, Grad)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2021 Pure Devices GmbH, Wuerzburg, Germany
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
Seq = set_EmptyField(Seq, 'nSaturationPulses', 1);
Seq = set_EmptyField(Seq, 'tFlip', 10e-3);
Seq = set_EmptyField(Seq, 'tFlipLog', 0);
if numel(Seq.tFlip) == 1
  Seq = set_EmptyField(Seq, 'tFlipStart', Seq.tFlip);
  Seq = set_EmptyField(Seq, 'tFlipEnd', Seq.tFlip*Seq.nFlips);
  % Create vector with times for the inversion pulse
  if Seq.tFlipLog
    Seq.tFlip = logspace(log10(Seq.tFlipStart), log10(Seq.tFlipEnd), Seq.nFlips);
  else
    Seq.tFlip = linspace(Seq.tFlipStart, Seq.tFlipEnd, Seq.nFlips);
  end
end
Seq = set_EmptyField(Seq, 'tRelax', Seq.LoopsBreak);
Seq = set_EmptyField(Seq, 'Recovery', 'Inversion');
Seq = set_EmptyField(Seq, 'tGradEC', 5e-3);
iDevice = Seq.AQSlice(1).iDevice;
Seq = set_EmptyField(Seq, 'tSaturationGrad', ...
  min(20e-3, Seq.tFlip(1)-Seq.tGradEC-HW.Grad(iDevice).tRamp*2-HW.Grad(iDevice).tEC*2-HW.tFlip90Def-0.5e-3));
Seq = set_EmptyField(Seq, 'AmpSaturationGrad', ...
  min(100e-3, min(HW.Grad(iDevice).MaxAmp(1:3))/5));


Seq.tFlip = Seq.tFlip([ones(1,Seq.nPreLoops),1:end]);  % repeat the first flip time for the pre-loops

%% Add pulses
switch Seq.Recovery
  case 'Inversion'
    Seq.iLaplace1D.Problem = 'Inversion';

    % 90 degrees slice excitation
    % Set the pulse shape for the initial inversion pulse:
    Seq.Slice(3).Pulse.Function = Seq.inversionPulse;
    Seq.Slice(3).Pulse.FlipAngleComposite = Seq.AQSlice(1).inversionFlipAngleComposite;
    Seq.Slice(3).Pulse.FlipPhaseComposite = Seq.AQSlice(1).inversionFlipPhaseComposite;
    % Make it an inversion pulse:
    Seq.Slice(3).Pulse.FlipAngle = 180;
    % Basically without slice selection
    Seq.Slice(3).Thickness = 1e12;
    % Set the center of the inversion pulse to "tFlip" seconds prior to the the
    % center of the excitation pulse. "Seq.Loop" is the loop counter. Thus, use
    % a higher tFlip every loop.
    Seq.Slice(3).CenterOfPulse = Seq.tEcho/2 - Seq.tFlip(Seq.Loop);
    % Add this pulse to the tRep with the excitation pulse(s):
    Seq.Slice(3).UseAtRepetitionTime = Seq.P90tReps;

    % Seq.Slice(3).UseCoordinate=Seq.AQSlice(1).SliceCoordinate;
    % Seq.Slice(3).GradTimeDelay=Seq.AQSlice(1).SliceTimeDelay;
    % Seq.Slice(3).Pulse.Phase=Seq.tRepTurboBlockPhaseExcitation(Seq.Slice(1).UseAtRepetitionTime)+Seq.AQSlice(1).excitationPhase;
    % Seq.Slice(3).Pulse.PhaseIncrement=Seq.AQSlice(1).excitationPhaseIncrement;
    % Seq.Slice(3).CenterOfRephase=Seq.Read(1).CenterOfRephase-Seq.tReadoutOffset-mean(Seq.Read(1).GradTimeDelayOffset);
    % Seq.Slice(3).GradRephaseLength=Seq.Read(1).GradRephaseLength;
    % Seq.Slice(3).MaxGradAmp=Seq.MaxGradAmpSlice;
    % Seq.Slice(3).distance=Seq.AQSlice(1).Center2OriginImage(1);
    % Seq.Slice(3).GradTimeIntegralRephaseOffset=Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset;
    % Seq.Slice(3).GradSign=Seq.AQSlice(1).SliceGradSign*Seq.AQSlice(1).GradSign;
    % Seq.Slice(3).GradDephaseSign=-Seq.AQSlice(1).SliceGradSign*Seq.AQSlice(1).GradSign;
    % Seq.Slice(3).GradRephaseSign=-Seq.AQSlice(1).SliceGradSign*Seq.AQSlice(1).GradSign;

    % generate slice parameters
    Seq = get_SliceParameter(Seq, HW);


    % 180 degrees spoiler encoding aligned with read rephase
    % Size of the spoiler in meter:
    Seq.Phase(10).sizePhase = Seq.InversionSpoilSize;
    % The spoiler does not actually encode the image in phase direction:
    Seq.Phase(10).nPhase = 1;
    % "Default" settings for a spoiler
    Seq.Phase(10).PhaseOS = 2;
    Seq.Phase(10).StepOrder = [1,1];
    % Length of the spoiler in seconds:
    Seq.Phase(10).GradDephaseLength = Seq.InversionSpoilLength;
    % Apply spoiler in read direction (1):
    Seq.Phase(10).UseCoordinate = 1;
    % Add the spoiler in the same tRep as the inversion pulse:
    Seq.Phase(10).UseAtRepetitionTime = Seq.Slice(3).UseAtRepetitionTime;
    % Sign of the spoiler gradient amplitude (arbitrary):
    Seq.Phase(10).GradDephaseSign = 1;
    Seq.Phase(10).CenterOfDephase = Seq.Slice(3).CenterOfPulse+Seq.Slice(3).GradLength/2+Seq.Phase(10).GradDephaseLength/2;
    % A rephase gradient before the inversion pulse is not necessary and it
    % won't be added to the sequence below. But add it for completeness anyway:
    Seq.Phase(10).GradRephaseLength = Seq.Phase(10).GradDephaseLength;
    Seq.Phase(10).GradRephaseSign = Seq.Phase(10).GradDephaseSign;
    Seq.Phase(10).CenterOfRephase = Seq.Slice(3).CenterOfPulse-Seq.Slice(3).GradLength/2-Seq.Phase(10).GradRephaseLength/2;

    % generate phase parameters
    Seq = get_PhaseParameter(Seq, HW);

    TX = add_TX(TX,Seq.Slice(3).TX);  % Add the inversion pulse
    Grad = add_Grad(Grad, Seq.Slice(3).Grad);  % Add the slice selection gradient
    Grad = add_Grad(Grad, Seq.Phase(10).GradDephase);  % Add the spoiler

  case 'Saturation'
    Seq.iLaplace1D.Problem = 'Saturation';
    Seq.tSaturat = HW.tFlip90Def*Seq.saturationPulse(HW,'Amp');
    % FIXME: Use newer calling convention for pulse shape functions.
    pulsePrepare = Seq.saturationPulse(HW, 0, 1/Seq.tSaturat*Seq.saturationPulse(HW,'Time'), pi*HW.tFlip90Def/(HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef), 51, Seq.tSaturat, HW.fLarmor, 0);

    % Gradients
    if Seq.tSaturat/2+HW.Grad(iDevice).tEC+Seq.tSaturationGrad+HW.Grad(iDevice).tRamp+HW.Grad(iDevice).tEC+Seq.tSaturat/2+50e-6 > Seq.tFlip(1)
      error('Increase Seq.tFlipStart or decrease Seq.tSaturationGrad');
    end
    temp = nan(size(Seq.tRep));
    temp(Seq.P90tReps)=1;
    for t = 1:HW.Grad(iDevice).n
      if t <= 3
        Seq.GradRecovery(t).Time = ...
          cumsum([Seq.tSaturat/2+HW.Grad(iDevice).tEC; ...
          HW.Grad(iDevice).tRamp; ...
          Seq.tSaturationGrad-HW.Grad(iDevice).tRamp*1; ...
          HW.Grad(iDevice).tRamp]) * temp - ...
          Seq.tFlip(Seq.Loop) + Seq.Slice(1).CenterOfPulse;
        Seq.GradRecovery(t).Amp = Seq.AmpSaturationGrad * [0;1;1;0;] * temp;
      else
        Seq.GradRecovery(t).Time = NaN;
        Seq.GradRecovery(t).Amp = NaN;
      end
    end

    if Seq.GradRecovery(1).Time(4,1) > Seq.Slice(1).CenterOfPulse-Seq.Slice(1).GradLength/2-Seq.tGradEC
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
      repmat(pulsePrepare.Start-Seq.tFlip(Seq.Loop)+Seq.Slice(1).CenterOfPulse, 1, numel(Seq.P90tReps));
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

