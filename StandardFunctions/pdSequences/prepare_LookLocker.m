function [HW, Seq, AQ, TX, Grad] = prepare_LookLocker(HW, Seq, AQ, TX, Grad)
%% Add inversion pulses into gradient echo sequence to measure T1* time
%
%     [HW, Seq, AQ, TX, Grad] = prepare_LookLocker(HW, Seq, AQ, TX, Grad)
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
%     preparationPulse
%         Pulse shape function of the 180 degrees preparation pulse.
%         (Default: @Pulse_Rect_Composite180)
%
%     preparationPulseProperties
%         Additional properties that should be passed to the pulse shape
%         function.
%
%     preparationSpoilSize
%         Size in meters of the spoiler after the preparation pulse.
%         (Default: 1e-3)
%
%     preparationSpoilLength
%         Duration in seconds of the spoiler after the preparation pulse.
%         (Default: 2e-3)
%
%     tPrepare
%         Time between the center of the preparation pulse and the center of the
%         next excitation pulse in seconds.
%         (Default: 1.5)
%
%     tRelax
%         Relaxation time in seconds before the next preparation pulse.
%         (Default: 1.5)
%
%     Recovery
%         Type of the preparation of spin system. The only supported type is
%         'Inversion'.
%         (Default: 'Inversion')
%
%     tShiftPreparation
%         Optionally, shift all preparation pulses by a fixed duration in
%         seconds. This value adds to the recovery time.
%         (Default: 0)
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
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% early return
if ~strcmp(Seq.LoopName{Seq.LoopNameCount}, 'normal')
  % only change normal loops
  return;
end

%% Default parameters
if isemptyfield(Seq, 'preparationPulse')
  Seq.preparationPulse = @Pulse_Rect_Composite180;
end
if isemptyfield(Seq, 'preparationPulseProperties')
  Seq.preparationPulseProperties = struct();
end
if isemptyfield(Seq, 'tPrepare')
  Seq.tPrepare = 10e-3;
end
if isemptyfield(Seq, 'tShiftPreparation')
  Seq.tShiftPreparation = 0;
end
if isemptyfield(Seq, 'tRelax')
  Seq.tRelax = 1.5;
end
if isemptyfield(Seq, 'Recovery')
  Seq.Recovery = 'Inversion';
end
% iDevice = Seq.AQSlice(1).iDevice;


%% Add pulses
switch Seq.Recovery
  case 'Inversion'
    Seq.iLaplace1D.Problem = 'Inversion';

    % default parameters specific to inversion recovery
    if isemptyfield(Seq, 'preparationSpoilSize')
      % size of spoiler after preparation pulse in meter
      Seq.preparationSpoilSize = 1e-3;
    end
    if isemptyfield(Seq, 'preparationSpoilLength')
      % length of spoiler after preparation pulse in seconds
      Seq.preparationSpoilLength = 2e-3;
    end

    % find tReps with first excitation for each line in k-space
    idxPrep = Seq.SteadyState_PreShots+1:Seq.AQSlice(1).nImagesSlice:numel(Seq.Slice(1).UseAtRepetitionTime);
    if Seq.SteadyState_PreShots >= Seq.AQSlice(1).nImagesSlice
      % insert preparation pulses in pre-shots
      idxPrep = [fliplr(Seq.SteadyState_PreShots+1-Seq.AQSlice(1).nImagesSlice:-Seq.AQSlice(1).nImagesSlice:1), ...
        idxPrep];
    end
    tRepkLines = Seq.Slice(1).UseAtRepetitionTime(idxPrep);

    % extend previous tRep to make room for preparation pulse
    if tRepkLines(1) == 1
      Seq.tRep(tRepkLines(2:end)-1) = Seq.tRep(tRepkLines(2:end)-1) + Seq.tRelax + Seq.tPrepare;
    else
      Seq.tRep(tRepkLines-1) = Seq.tRep(tRepkLines-1) + Seq.tRelax + Seq.tPrepare;
    end

    tSlice = numel(Seq.Slice);
    % 90 degrees slice excitation
    % Set the pulse shape for the initial inversion pulse:
    Seq.Slice(tSlice+1).Pulse.Function = Seq.preparationPulse;
    % Make it an inversion pulse:
    Seq.Slice(tSlice+1).Pulse.FlipAngle = 180;
    % use phase of next excitation pulse as reference
    Seq.Slice(tSlice+1).Pulse.Phase = Seq.Slice(1).Pulse.Phase(idxPrep) + 90;
    % without slice selection
    % FIXME: Should we support slice selective inversion pulses?
    Seq.Slice(tSlice+1).Thickness = Inf;
    % Set the center of the inversion pulse to "tFlip" seconds prior to the the
    % center of the excitation pulse. "Seq.Loop" is the loop counter. Thus, use
    % a higher tFlip every loop.
    Seq.Slice(tSlice+1).CenterOfPulse = Seq.Slice(1).CenterOfPulse - Seq.tPrepare - Seq.tShiftPreparation;
    % Add this pulse to the tRep with the excitation pulse(s):
    Seq.Slice(tSlice+1).UseAtRepetitionTime = tRepkLines;
    % set optional properties for pulse shape functions
    if ~isemptyfield(Seq.AQSlice(1), 'inversionFlipAngleComposite')
      Seq.Slice(tSlice+1).Pulse.FlipAngleComposite = Seq.AQSlice(1).inversionFlipAngleComposite;
    end
    if ~isemptyfield(Seq.AQSlice(1), 'inversionFlipPhaseComposite')
      Seq.Slice(tSlice+1).Pulse.FlipPhaseComposite = Seq.AQSlice(1).inversionFlipPhaseComposite;
    end
    preparationPulseProperties = fieldnames(Seq.preparationPulseProperties);
    for iProp = 1:numel(preparationPulseProperties)
      % copy additional properties for pulse shape function
      Seq.Slice(tSlice+1).Pulse.(preparationPulseProperties{iProp}) = ...
        Seq.preparationPulseProperties.(preparationPulseProperties{iProp});
    end

    % Seq.Slice(tSlice+1).UseCoordinate=Seq.AQSlice(1).SliceCoordinate;
    % Seq.Slice(tSlice+1).GradTimeDelay=Seq.AQSlice(1).SliceTimeDelay;
    % Seq.Slice(tSlice+1).Pulse.PhaseIncrement = Seq.AQSlice(1).excitationPhaseIncrement;
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


    % spoiler of 180 degrees preparation pulse
    tPhase = numel(Seq.Phase);
    % size of spoiler in meter:
    Seq.Phase(tPhase+1).sizePhase = Seq.preparationSpoilSize;
    % Note: spoiler does not actually encode image in phase direction
    Seq.Phase(tPhase+1).nPhase = 1;
    % "default" settings for a spoiler
    Seq.Phase(tPhase+1).PhaseOS = 2;
    Seq.Phase(tPhase+1).StepOrder = [1, 1];
    % length of spoiler in seconds:
    Seq.Phase(tPhase+1).GradDephaseLength = Seq.preparationSpoilLength;
    % apply spoiler in slice direction
    Seq.Phase(tPhase+1).UseCoordinate = Seq.AQSlice(1).SliceCoordinate;
    % add spoiler in same tRep as the inversion pulse:
    Seq.Phase(tPhase+1).UseAtRepetitionTime = Seq.Slice(tSlice+1).UseAtRepetitionTime;
    % sign of the spoiler gradient amplitude (arbitrary):
    Seq.Phase(tPhase+1).GradDephaseSign = 1;
    Seq.Phase(tPhase+1).CenterOfDephase = Seq.Slice(tSlice+1).CenterOfPulse + Seq.Slice(tSlice+1).GradLength/2 + Seq.Phase(tPhase+1).GradDephaseLength/2;
    % A rephase gradient before the inversion pulse is not necessary and it
    % won't be added to the sequence below. But add it for completeness anyway:
    Seq.Phase(tPhase+1).GradRephaseLength = Seq.Phase(tPhase+1).GradDephaseLength;
    Seq.Phase(tPhase+1).GradRephaseSign = Seq.Phase(tPhase+1).GradDephaseSign;
    Seq.Phase(tPhase+1).CenterOfRephase = Seq.Slice(tSlice+1).CenterOfPulse - Seq.Slice(tSlice+1).GradLength/2 - Seq.Phase(tPhase+1).GradRephaseLength/2;

    % generate phase parameters
    Seq = get_PhaseParameter(Seq, HW);

    TX = add_TX(TX, Seq.Slice(tSlice+1).TX);  % add inversion pulse (for preparation)
    % Grad = add_Grad(Grad, Seq.Slice(tSlice+1).Grad);  % add slice selection gradient
    Grad = add_Grad(Grad, Seq.Phase(tPhase+1).GradDephase);  % add spoiler (after preparation pulse)

  otherwise
    error('Seq.Recovery must be "Inversion"!');
    % FIXME: Do other types of preparation even make sense?

end

end
