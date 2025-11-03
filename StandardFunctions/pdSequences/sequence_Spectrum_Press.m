function [SeqLoop, mySave] = sequence_Spectrum_Press(HW, Seq, AQ, TX, Grad, mySave)
%% Acquire a signal using a PRESS sequence
%
%   [SeqLoop, mySave] = sequence_Spectrum_Press(HW, Seq, AQ, TX, Grad, mySave)
%
% A PRESS (Point RESolved Sprectroscopy) sequence is a multi spin echo sequence
% which acquires the second echo. The excitation and refocusing pulses are slice
% selective in orientations that are orthogonal to each other.
%
% This function is similar to "sequence_Spin_Echo". Many settings of that
% function also apply to this function. Only the most notable differences are
% documented here.
%
% The acquisition starts at the center of the second echo.
%
% The correction loops of "sequence_Spin_Echo" do not apply.
%
% The default parameters are set for water.
%
% INPUT:
%
%   HW
%           HW object or structure (see LoadSystem).
%
%   Seq
%           Structure containing the settings for the sequence. All settings are
%           optional and default values are used if they are omitted. Among
%           others the following fields can be set:
%
%     Spectro
%           Structure with spectroscopy settings. All settings are optional and
%           default values are used if they are omitted. Among others the
%           following fields can be set:
%
%       fOffsetTimeStart
%             Especially important when averaging, the frequency of the acquired
%             signal is determined to compensate the B0 change potentially
%             caused by temperature drift of the magnet. The value determined
%             for each averaging step (Seq.Loop) is returned in
%             "SeqLoop.loopdata.fOffset". This value defines the offset in
%             seconds with respect to the echo center that corresponds to the
%             start of the interval used for frequency determination.
%             (Default: 0.03)
%
%       fOffsetTimeStop
%             Corresponds to "Seq.Spectro.fOffsetTimeStart" but defines the
%             offset in seconds with respect to the echo center that corresponds
%             to the end of the interval used for frequency determination.
%             (Default: 0.09)
%
%       PhaseOffsetTimeStart
%             The temperature drift might also cause a phase shift of the
%             acquired signal. The average phase of the signal is determined for
%             each averaging step (Seq.Loop) and returned in
%             "SeqLoop.loopdata.phaseOffset". This value defines the offset in
%             seconds with respect to the echo center that corresponds to the
%             start of the interval used for signal phase determination.
%             (Default: 0.03)
%
%       PhaseOffsetTimeStop
%             Corresponds to "Seq.Spectro.PhaseOffsetTimeStart" but defines the
%             offset in seconds with respect to the echo center that corresponds
%             to the end of the interval used for signal phase determination.
%             (Default: 0.04)
%
%     AQSlice
%           scalar structure defining the orientation and other properties of
%           the acquired slice. Among others the following fields can be used:
%
%       thickness
%             See the documentation for "sequence_Spin_Echo".
%
%       thicknessInversion
%             See the documentation for "sequence_Spin_Echo".
%
%       thicknessInversion2
%             Like "Seq.AQSlice.thicknessInversion" but applies to the second
%             inversion pulse.
%
%       excitationPulse
%             See the documentation for "sequence_Spin_Echo".
%
%       inversionPulse
%             See the documentation for "sequence_Spin_Echo".
%
%       inversionPulse2
%             Like "Seq.AQSlice.inversionPulse" but applies to the second
%             inversion pulse.
%
%   AQ, TX, Grad, mySave
%         See documentation for "sequence_Spin_Echo".
%
%
% OUTPUT:
%
%   Similar to "sequence_Spin_Echo".
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2021 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------


% if nargin<4
%   if exist(HW.Spectro.ReferenceFidPath, 'file')
%     if HW.Spectro.useSliceSelect  % use slice gradient
%       Seq.useSliceSelect = 1;
%       SliceSelect = HW.Spectro.SliceSelect;
%       Seq.SlicePulse = HW.Spectro.SlicePulse;
%     else
%       Seq.useSliceSelect = 0;
%       SliceSelect = [];
%       Seq.FlipPulse = HW.Spectro.FlipPulse;
%     end
%   end
% end

if isemptyfield(Seq, 'Loops'), Seq.Loops = 1; end  % 1,2,3....
if isemptyfield(Seq, 'PreLoops'), Seq.PreLoops = 1; end  % 1,2,3....
if isemptyfield(Seq, 'plot'), Seq.plot = 0; end  % 1,2,3....
if isemptyfield(Seq, 'fSample'), Seq.fSample = HW.FindShim.fSample; end  % HW.RX.fSample = 125e6 Hz
if isemptyfield(Seq, 'fSampleFID'), Seq.fSampleFID = Seq.fSample; end  % HW.RX.fSample = 125e6 Hz
if isemptyfield(Seq, 'LoopsRepetitionTime'), Seq.LoopsRepetitionTime = HW.FindFrequencyPause; end
if isemptyfield(Seq, 'Find_Frequency_interval'), Seq.Find_Frequency_interval = 0; end  % Find Frequency interval

if isemptyfield(Seq, 'Spectro'), Seq.Spectro = []; end
if isemptyfield(Seq.Spectro, 'fOffsetTimeStart'), Seq.Spectro.fOffsetTimeStart = 0.03; end
if isemptyfield(Seq.Spectro, 'fOffsetTimeStop'), Seq.Spectro.fOffsetTimeStop = 0.09; end
if isemptyfield(Seq.Spectro, 'PhaseOffsetTimeStart'), Seq.Spectro.PhaseOffsetTimeStart = 0.03; end
if isemptyfield(Seq.Spectro, 'PhaseOffsetTimeStop'), Seq.Spectro.PhaseOffsetTimeStop = 0.04; end

if nargin < 3
  AQ.Start = [];        % matrix containing the starting times (center of the first Sample – 0.5/AQ.fSample) of each acquisition window for each TR; up to 510 AQ windows are possible per TR
  AQ.nSamples = [];     % matrix containing the number of samples to be acquired for each AQ window for each TR
  AQ.fSample = [];      % matrix containing the sampling frequency of each acquisition window for each TR
  AQ.Frequency = [];    % matrix containing the mixing frequency of each acquisition windows for each TR
  AQ.Phase = [];        % matrix containing the phase offset for each acquisition window
  AQ.ResetPhases = [];  % if set to 1 the TX and RX phase is matched at the beginning of the TR
  AQ.Gain = [];         % relative gain of the acquisition path (1/AQ.Gain = max input Amplitude). Use HW.RX.GainDef for the best noise figure.
  AQ.Repeat = [];       % if there is a ‘1’ in the array, the settings for the window of the prior TR are used. (reduces data traffic between the PC and the MRI device)
end

if nargin < 4
  TX.Channel = HW.TX.ChannelDef;  % use transmitting channel 1 or 2
  TX.Start = [];        % matrix containing the starting time of the RF pulses for each TR column by column; up to 510 RF pulses are possible per TR (column)
  TX.Frequency = [];    % matrix containing transmitting frequency for each pulse of each TR
  TX.Duration = [];     % matrix containing the duration of the RF pulses for each TR column by column; up to 510 RF pulses are possible per TR
  TX.Amplitude = [];    % matrix containing the amplitude for each pulse of each TR
  TX.Phase = [];        % matrix containing the phase offset for each pulse of each TR

  TX.BlankOffset = [];  % array of times to define, when the blanking signal is applied prior to the RF pulse; can be adjusted every TR
  TX.BlankPostset = []; % array of times to define how long blanking stays active after a RF pulse; can be adjusted every TR
  TX.Repeat = [];       % if there is a ‘1’ in the array, the pulses of the prior TR are used. (reduces data traffic between the PC and the MRI device)
end

iDevice = 1;  % FIXME: Support multiple MMRT devices
Channel = 1;  % FIXME: Add support for multiple AQ channels?

if nargin < 5
  [Grad(1:HW.Grad(iDevice).n).Time] = deal([]);   % time the values in Grad(t).Amp with the corresponding index are set
  [Grad(1:HW.Grad(iDevice).n).Amp] = deal([]);    % amplitude of the gradient at the time Grad(t).Time (linearly interpolated)
  [Grad(1:HW.Grad(iDevice).n).Shim] = deal([]);   % additional shim; magnet shim is already considered in HW.MagnetShim. Caution: Using high values over a long time will damage the gradient coils and amplifiers!
  [Grad(1:HW.Grad(iDevice).n).Repeat] = deal([]); % if there is a one in the array the gradients of the prior TR are used. (reduces data traffic between the PC and the MRI device)
end

if nargin >= 6
  HW.Grad(iDevice).HoldShim = 1;
  HW.FindFrequencyPause = Seq.LoopsRepetitionTime;
  % update approximate magnet frequency if 1000s have passed since last frequency update
  [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 1000);
  % repeat frequency sweep with improved accuracy
  [HW, mySave] = Find_Frequency_Sweep(HW, mySave, Seq.Find_Frequency_interval, ...
    0, 1, HW.tFlip90Def*HW.FindFrequencySweep.tPulseDivider, 1, 1024);
else
  mySave = [];
end


if (Seq.LoopsRepetitionTime <= 30) && ((Seq.Loops+Seq.PreLoops) > 1)
  Seq.TimeToNextSequence = Seq.LoopsRepetitionTime;
end

tc = 0;
PreLoops = ones(1, Seq.PreLoops);

for t = 1:Seq.PreLoops
  tc = tc+1;
  Seq.LoopName{tc} = 'normal PreLoops';
end
for t = 1:Seq.Loops
  tc = tc+1;
  Seq.LoopName{tc} = 'normal';
end


if isemptyfield(Seq, 'tFID'), Seq.tFID = 1; end  % duration of FID acquisition in spectroscopy experiment
if isemptyfield(Seq, 'tEcho'), Seq.tEcho = 2e-3; end
if isemptyfield(Seq, 'tRep'), Seq.tRep = Seq.tEcho*3 + Seq.tFID + 1e-3; end
if isemptyfield(Seq, 'nEchos'), Seq.nEchos = 0; end
if isemptyfield(Seq, 'Reinitialize'), Seq.Reinitialize = 1; end
if ~isfield(Seq, 'AQSlice'), Seq.AQSlice = []; end
if isemptyfield(Seq.AQSlice, 'excitationPulse'), Seq.AQSlice(1).excitationPulse = @Pulse_RaisedCos; end
if isemptyfield(Seq.AQSlice, 'inversionPulse'), Seq.AQSlice(1).inversionPulse = Seq.AQSlice(1).excitationPulse; end
if isemptyfield(Seq.AQSlice, 'inversionPulse2'), Seq.AQSlice(1).inversionPulse2 = Seq.AQSlice(1).inversionPulse; end
if isemptyfield(Seq.AQSlice, 'excitationFlipAngle'), Seq.AQSlice(1).excitationFlipAngle = 90; end
if isemptyfield(Seq.AQSlice, 'inversionFlipAngle'), Seq.AQSlice(1).inversionFlipAngle = 180; end
if ~isfield(Seq.AQSlice, 'excitationFlipAngleComposite'), Seq.AQSlice(1).excitationFlipAngleComposite = []; end
if ~isfield(Seq.AQSlice, 'excitationFlipPhaseComposite'), Seq.AQSlice(1).excitationFlipPhaseComposite = []; end
if ~isfield(Seq.AQSlice, 'inversionFlipAngleComposite'), Seq.AQSlice(1).inversionFlipAngleComposite = []; end
if ~isfield(Seq.AQSlice, 'inversionFlipPhaseComposite'), Seq.AQSlice(1).inversionFlipPhaseComposite = []; end
if isemptyfield(Seq.AQSlice, 'inversion2FlipAngleComposite'), Seq.AQSlice(1).inversion2FlipAngleComposite = Seq.AQSlice(1).inversionFlipAngleComposite; end
if isemptyfield(Seq.AQSlice, 'inversion2FlipPhaseComposite'), Seq.AQSlice(1).inversion2FlipPhaseComposite = Seq.AQSlice(1).inversionFlipAngleComposite; end
if isemptyfield(Seq.AQSlice, 'excitationPhase'), Seq.AQSlice(1).excitationPhase = 0; end
if isemptyfield(Seq.AQSlice, 'inversionPhase'), Seq.AQSlice(1).inversionPhase = 90; end
if isemptyfield(Seq.AQSlice, 'excitationPhaseIncrement'), Seq.AQSlice(1).excitationPhaseIncrement = 0; end
if isemptyfield(Seq.AQSlice, 'inversionPhaseIncrement'), Seq.AQSlice(1).inversionPhaseIncrement = 0; end
if isemptyfield(Seq.AQSlice, 'thickness'), Seq.AQSlice(1).thickness = 5e-3; end
if isemptyfield(Seq.AQSlice, 'thicknessInversion'), Seq.AQSlice(1).thicknessInversion = Seq.AQSlice(1).thickness; end
if isemptyfield(Seq.AQSlice, 'thicknessInversion2'), Seq.AQSlice(1).thicknessInversion2 = Seq.AQSlice(1).thicknessInversion; end
if isemptyfield(Seq.AQSlice, 'MaxGradAmpSlice'), Seq.AQSlice(1).MaxGradAmpSlice = HW.Grad.MaxAmpSlice; end
if isemptyfield(Seq.AQSlice, 'MaxGradAmpInversion'), Seq.AQSlice(1).MaxGradAmpInversion = Seq.AQSlice(1).MaxGradAmpSlice; end
if isemptyfield(Seq.AQSlice, 'Center2OriginImage'), Seq.AQSlice(1).Center2OriginImage = [0, 0, 0]; end
if isemptyfield(Seq.AQSlice, 'sizePhaseSpoil'), Seq.AQSlice(1).sizePhaseSpoil = [Inf, Inf, Inf, Inf, Inf]; end
if isemptyfield(Seq.AQSlice, 'ReadOS'), Seq.AQSlice(1).ReadOS = 1; end
if isemptyfield(Seq.AQSlice, 'nRead'), Seq.AQSlice(1).nRead = round(Seq.tFID*Seq.fSampleFID); end
if isemptyfield(Seq.AQSlice, 'sizeRead'), Seq.AQSlice(1).sizeRead = Inf; end
if isemptyfield(Seq.AQSlice, 'alfa'), Seq.AQSlice(1).alfa = 0; end
if isemptyfield(Seq.AQSlice, 'phi'), Seq.AQSlice(1).phi = 0; end
if isemptyfield(Seq.AQSlice, 'theta'), Seq.AQSlice(1).theta = 0; end
if isemptyfield(Seq.AQSlice, 'tCrusherSliceEqual'), Seq.AQSlice(1).tCrusherSliceEqual = 0; end
if isemptyfield(Seq.AQSlice, 'tCrusher'), Seq.AQSlice(1).tCrusher = 1e-3; end
if isemptyfield(Seq.AQSlice, 'tCrusherEC'), Seq.AQSlice(1).tCrusherEC = HW.Grad(iDevice).tEC; end
if isemptyfield(Seq.AQSlice, 'tCrusherAdd'), Seq.AQSlice(1).tCrusherAdd = Seq.tEcho/5; end
if isemptyfield(Seq.AQSlice, 'StartWithKSpaceCenter'), Seq.AQSlice(1).StartWithKSpaceCenter = 1; end
if isemptyfield(Seq.AQSlice, 'ReadCoordinate'), Seq.AQSlice(1).ReadCoordinate = 1; end
if isemptyfield(Seq, 'LoopsBreakExactly'), Seq.LoopsBreakExactly = 0; end
if ~isfield(Seq, 'LoopsBreak'), Seq.LoopsBreak = []; end


%% initial values for all loops
init.AQ = AQ;
init.TX = TX;
init.Seq = Seq;
init.Grad = Grad;
LoopCount = 0;

%% loops
for Loop = [PreLoops, 1:Seq.Loops]
  LoopCount = LoopCount+1;
  AQ = init.AQ;
  TX = init.TX;
  Grad = init.Grad;
  Seq = init.Seq;
  Seq.Loop = Loop;
  Seq.LoopNameCount = LoopCount;


  if Seq.Loops+Seq.PreLoops>1
    disp(['loops to run ' num2str(Seq.Loops-Loop+1) ' (' num2str((Seq.Loops-Loop+1)*Seq.LoopsRepetitionTime) ' s)' ]);
  end


  % default parameters
  Seq = set_EmptyField(Seq, 'PreProcessSequence', 1);
  Seq = set_EmptyField(Seq, 'StartSequence', 1);
  Seq = set_EmptyField(Seq, 'PollPPGfast', 1);
  Seq = set_EmptyField(Seq, 'GetRawData', 1);
  Seq = set_EmptyField(Seq, 'PostProcessSequence', 1);


  % Readout at tEcho
  Seq.Read(1).HzPerPixelMin = 1/Seq.tFID;
  Seq.Read(1).CenterOfReadout = 0;
  Seq.Read(1).UseCoordinate = Seq.AQSlice(1).ReadCoordinate;
  Seq.Read(1).GradTimeIntegralOffset = 0;
  Seq.Read(1).UseAtRepetitionTime = 1;
  Seq.Read(1).Phase = 0;
  Seq.Read(1).PhaseIncrement = 0;
  Seq.Read(1).distance = 0;
  Seq.Read(1).GradSign = 1;
  Seq.Read(1).useAQSlice = 1;
  Seq.Read(1).StartWithKSpaceCenter = Seq.AQSlice(1).StartWithKSpaceCenter;


  % 90 degrees slice excitation
  Seq.Slice(1).Pulse.Function = Seq.AQSlice(1).excitationPulse;
  Seq.Slice(1).Pulse.FlipAngleComposite = Seq.AQSlice(1).excitationFlipAngleComposite;
  Seq.Slice(1).Pulse.FlipPhaseComposite = Seq.AQSlice(1).excitationFlipPhaseComposite;
  Seq.Slice(1).Pulse.FlipAngle = Seq.AQSlice(1).excitationFlipAngle;
  Seq.Slice(1).Thickness = Seq.AQSlice(1).thickness;
  Seq.Slice(1).CenterOfPulse = -4*Seq.tEcho/2;
  Seq.Slice(1).UseCoordinate = 1;
  Seq.Slice(1).UseAtRepetitionTime = 1;
  Seq.Slice(1).Pulse.Phase = Seq.AQSlice(1).excitationPhase;
  Seq.Slice(1).Pulse.PhaseIncrement = Seq.AQSlice(1).excitationPhaseIncrement;
  % Seq.Slice(1).CenterOfRephase = Seq.Slice(1).CenterOfPulse+Seq.tEcho/4;
  Seq.Slice(1).GradRephaseLength = Seq.AQSlice(1).tCrusherAdd;
  Seq.Slice(1).MaxGradAmp = Seq.AQSlice(1).MaxGradAmpSlice;
  Seq.Slice(1).distance = Seq.AQSlice(1).Center2OriginImage(1);
  % Seq.Slice(1).GradTimeIntegralRephaseOffset = 0;
  Seq.Slice(1).GradSign = 1;
  Seq.Slice(1).GradDephaseSign = -1;
  Seq.Slice(1).GradRephaseSign = -1;
  Seq.Slice(1).useAQSlice = 1;


  % 180 degrees (slice) inversion
  Seq.Slice(2).Pulse.Function = Seq.AQSlice(1).inversionPulse2;
  Seq.Slice(2).Pulse.FlipAngleComposite = Seq.AQSlice(1).inversionFlipAngleComposite;
  Seq.Slice(2).Pulse.FlipPhaseComposite = Seq.AQSlice(1).inversionFlipPhaseComposite;
  Seq.Slice(2).Pulse.FlipAngle = Seq.AQSlice(1).inversionFlipAngle;
  Seq.Slice(2).Thickness = Seq.AQSlice(1).thicknessInversion;
  Seq.Slice(2).CenterOfPulse = -3*Seq.tEcho/2;
  Seq.Slice(2).UseCoordinate = 2;
  Seq.Slice(2).UseAtRepetitionTime = 1;
  Seq.Slice(2).Pulse.Phase = Seq.AQSlice(1).inversionPhase;
  Seq.Slice(2).Pulse.PhaseIncrement = Seq.AQSlice(1).inversionPhaseIncrement;
  % Seq.Slice(2).CenterOfRephase = Seq.Slice(2).CenterOfPulse+Seq.tEcho/4;
  Seq.Slice(2).GradRephaseLength = Seq.AQSlice(1).tCrusherAdd;
  % Seq.Slice(2).CenterOfDephase = Seq.Slice(2).CenterOfPulse-Seq.tEcho/4;
  Seq.Slice(2).GradDephaseLength = Seq.AQSlice(1).tCrusherAdd;
  Seq.Slice(2).MaxGradAmp = Seq.AQSlice(1).MaxGradAmpInversion;
  Seq.Slice(2).distance = Seq.AQSlice(1).Center2OriginImage(2);
  Seq.Slice(2).GradSign = 1;
  Seq.Slice(2).GradDephaseSign = -1;
  Seq.Slice(2).GradRephaseSign = -1;
  Seq.Slice(2).useAQSlice = 1;

  % 180 degrees (slice) inversion2
  Seq.Slice(3).Pulse.Function = Seq.AQSlice(1).inversionPulse;
  Seq.Slice(3).Pulse.FlipAngleComposite = Seq.AQSlice(1).inversionFlipAngleComposite;
  Seq.Slice(3).Pulse.FlipPhaseComposite = Seq.AQSlice(1).inversionFlipPhaseComposite;
  Seq.Slice(3).Pulse.FlipAngle = Seq.AQSlice(1).inversionFlipAngle;
  Seq.Slice(3).Thickness = Seq.AQSlice(1).thicknessInversion2;
  Seq.Slice(3).CenterOfPulse = -1*Seq.tEcho/2;
  Seq.Slice(3).UseCoordinate = 3;
  Seq.Slice(3).UseAtRepetitionTime = 1;
  Seq.Slice(3).Pulse.Phase = Seq.AQSlice(1).inversionPhase;
  Seq.Slice(3).Pulse.PhaseIncrement = Seq.AQSlice(1).inversionPhaseIncrement;
  % Seq.Slice(3).CenterOfRephase = Seq.Slice(3).CenterOfPulse+Seq.tEcho/4;
  Seq.Slice(3).GradRephaseLength = Seq.AQSlice(1).tCrusherAdd;
  % Seq.Slice(3).CenterOfDephase = Seq.Slice(3).CenterOfPulse-Seq.tEcho/4;
  Seq.Slice(3).GradDephaseLength = Seq.AQSlice(1).tCrusherAdd;
  Seq.Slice(3).MaxGradAmp = Seq.AQSlice(1).MaxGradAmpInversion;
  Seq.Slice(3).distance = Seq.AQSlice(1).Center2OriginImage(3);
  Seq.Slice(3).GradSign = 1;
  Seq.Slice(3).GradDephaseSign = -1;
  Seq.Slice(3).GradRephaseSign = -1;
  Seq.Slice(3).useAQSlice = 1;

  % generate readout timings
  Seq = get_SliceParameter(Seq, HW);

  % Seq.Slice(1).CenterOfRephase=Seq.Slice(2).CenterOfDephase;
  Seq.Slice(2).tEC = max(Seq.AQSlice(1).tCrusherSliceEqual/2-Seq.Slice(2).GradLength/2+Seq.Slice(2).tRamp, Seq.AQSlice(1).tCrusherEC);
  Seq.Slice(3).tEC = max(Seq.AQSlice(1).tCrusherSliceEqual/2-Seq.Slice(3).GradLength/2+Seq.Slice(3).tRamp, Seq.AQSlice(1).tCrusherEC);
  Seq.Slice(1).CenterOfRephase = [];
  Seq.Slice(1).CenterOfDephase = [];
  Seq.Slice(2).CenterOfRephase = [];
  Seq.Slice(2).CenterOfDephase = [];
  Seq.Slice(3).CenterOfRephase = [];
  Seq.Slice(3).CenterOfDephase = [];

  % generate readout timings
  Seq = get_SliceParameter(Seq, HW);
  Seq.Slice(1).CenterOfRephase = Seq.Slice(2).CenterOfDephase;
  Seq.Slice(1).CenterOfDephase = Seq.Slice(3).CenterOfRephase;


  % 90 degrees spoiler encoding aligned with read Seq.Slice(1).GradRephase and Seq.Slice(3).GradRephase
  Seq.Phase(1).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(1);
  Seq.Phase(1).nPhase = 1;
  Seq.Phase(1).PhaseOS = 2;
  Seq.Phase(1).StepOrder = [1, 1];
  Seq.Phase(1).CenterOfDephase = Seq.Slice(1).CenterOfRephase;
  Seq.Phase(1).CenterOfRephase = Seq.Slice(3).CenterOfRephase;
  Seq.Phase(1).GradDephaseLength = Seq.Slice(1).GradRephaseLength;
  Seq.Phase(1).GradRephaseLength = Seq.Slice(3).GradRephaseLength;
  Seq.Phase(1).UseCoordinate = 1;
  Seq.Phase(1).UseAtRepetitionTime = 1;
  Seq.Phase(1).GradDephaseSign = -1;
  Seq.Phase(1).GradRephaseSign = 1;

  % 180 degrees spoiler encoding aligned with read rephase
  Seq.Phase(2).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(2);
  Seq.Phase(2).nPhase = 1;
  Seq.Phase(2).PhaseOS = 2;
  Seq.Phase(2).StepOrder = [1, 1];
  Seq.Phase(2).CenterOfDephase = Seq.Slice(2).CenterOfDephase;
  Seq.Phase(2).CenterOfRephase = Seq.Slice(2).CenterOfRephase;
  Seq.Phase(2).GradDephaseLength = Seq.Slice(2).GradDephaseLength;
  Seq.Phase(2).GradRephaseLength = Seq.Slice(2).GradRephaseLength;
  Seq.Phase(2).UseCoordinate = 2;
  Seq.Phase(2).UseAtRepetitionTime = 1;
  Seq.Phase(2).GradDephaseSign = -1;
  Seq.Phase(2).GradRephaseSign = -1;

  % 180 degrees spoiler encoding aligned with read rephase
  Seq.Phase(3).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(3);
  Seq.Phase(3).nPhase = 1;
  Seq.Phase(3).PhaseOS = 2;
  Seq.Phase(3).StepOrder = [1, 1];
  Seq.Phase(3).CenterOfDephase = Seq.Slice(3).CenterOfDephase;
  Seq.Phase(3).CenterOfRephase = Seq.Slice(3).CenterOfRephase;
  Seq.Phase(3).GradDephaseLength = Seq.Slice(3).GradDephaseLength;
  Seq.Phase(3).GradRephaseLength = Seq.Slice(3).GradRephaseLength;
  Seq.Phase(3).UseCoordinate = 3;
  Seq.Phase(3).UseAtRepetitionTime = 1;
  Seq.Phase(3).GradDephaseSign = -1;
  Seq.Phase(3).GradRephaseSign = -1;


  % generate readout timings
  Seq = get_ReadParameter(Seq, HW);
  % generate slice parameters
  Seq = get_SliceParameter(Seq, HW);
  % generate phase parameters
  Seq = get_PhaseParameter(Seq, HW);
  % Gradients of image encoding

  % AQ
  AQ = add_AQ(AQ, Seq.Read(1).AQ);

  % TX
  for t = 1:3
    TX = add_TX(TX, Seq.Slice(t).TX);
  end
  Grad = add_Grad(Grad, Seq.Read(1).GradDephase);
  Grad = add_Grad(Grad, Seq.Read(1).Grad);
  Grad = add_Grad(Grad, Seq.Read(1).GradRephase);


  % Grad = add_Grad(Grad, Seq.Slice(1).GradDephase);
  Grad = add_Grad(Grad, Seq.Slice(1).Grad);
  Grad = add_Grad(Grad, Seq.Slice(1).GradRephase);

  Grad = add_Grad(Grad, Seq.Slice(2).Grad);
  Grad = add_Grad(Grad, Seq.Slice(3).Grad);


  Grad = add_Grad(Grad, Seq.Slice(2).GradDephase);
  Grad = add_Grad(Grad, Seq.Slice(3).GradDephase);
  Grad = add_Grad(Grad, Seq.Slice(2).GradRephase);
  Grad = add_Grad(Grad, Seq.Slice(3).GradRephase);

  if Seq.AQSlice(1).sizePhaseSpoil(1) < 1
    Grad = add_Grad(Grad, Seq.Phase(1).GradDephase);
    Grad = add_Grad(Grad, Seq.Phase(1).GradRephase);
  end

  if Seq.AQSlice(1).sizePhaseSpoil(2) < 1
    Grad = add_Grad(Grad, Seq.Phase(2).GradRephase);
    Grad = add_Grad(Grad, Seq.Phase(2).GradDephase);
  end

  if Seq.AQSlice(1).sizePhaseSpoil(3) < 1
    Grad = add_Grad(Grad, Seq.Phase(3).GradDephase);
    Grad = add_Grad(Grad, Seq.Phase(3).GradRephase);
  end

  if Seq.Reinitialize
    Seq.Reinitialize = double((Seq.LoopNameCount==1) || ~Seq.LoopsBreakExactly);
  end

  if ~any([Seq.PreProcessSequence, ...
           Seq.StartSequence, ...
           Seq.PollPPGfast, ...
           Seq.GetRawData, ...
           Seq.PostProcessSequence])
    SeqLoop = Seq;
    if isobject(HW)
      SeqLoop.HW = HW.ToStruct();
    else
      SeqLoop.HW = HW;
    end
    SeqLoop.AQ = AQ;
    SeqLoop.TX = TX;
    SeqLoop.Grad = Grad;
    break
  end

  if Seq.PreProcessSequence
    tStartSequence = Seq.StartSequence;
    Seq.StartSequence = 0;
    Seq.PollPPGfast = 0;
    Seq.GetRawData = 0;
    Seq.PostProcessSequence = 0;

    if tStartSequence
      [~, SeqOut] = set_sequence(HW, Seq, AQ, TX, Grad);  % PreProcessSequence sequence
    else
      [~, SeqLoop] = set_sequence(HW, Seq, AQ, TX, Grad); % PreProcessSequence sequence only
      return;
    end
  end

  if ~isempty(SeqOut.LoopsRepetitionTime) && isempty(SeqOut.LoopsBreak)
    init.Seq.LoopsBreak = SeqOut.LoopsRepetitionTime-SeqOut.SequenceTime;
    SeqOut.LoopsBreak = init.Seq.LoopsBreak;
  end

  if SeqOut.LoopsBreakExactly
    if SeqOut.LoopNameCount == 1
      SeqOut.Reinitialize = 1;
      if numel(SeqOut.LoopName) > 1
        SeqOut.TimeToNextSequence = SeqOut.LoopsBreak+SeqOut.tOffset(1);
        SeqOut.TimeFromLastSequence = [];
      else
        SeqOut.TimeToNextSequence = [];
        SeqOut.TimeFromLastSequence = [];
      end
    elseif (SeqOut.LoopNameCount>=2) && (SeqOut.LoopNameCount~=numel(SeqOut.LoopName))
      SeqOut.Reinitialize = 0;
      SeqOut.TimeFromLastSequence = SeqOut.LoopsBreak+SeqOut.tOffset(1);
      SeqOut.TimeToNextSequence = SeqOut.LoopsBreak+SeqOut.tOffset(1);
    elseif SeqOut.LoopNameCount == numel(SeqOut.LoopName)
      SeqOut.Reinitialize = 0;
      SeqOut.TimeFromLastSequence = SeqOut.LoopsBreak+SeqOut.tOffset(1);
      SeqOut.TimeToNextSequence = [];
    end
  end

  if isempty(SeqOut.LoopsBreak), tb=1; else tb=SeqOut.LoopsBreak; end
  if strcmp(SeqOut.LoopName{SeqOut.LoopNameCount}, 'normal') && (SeqOut.Loops > 1)
    if isempty(SeqOut.LoopsRepetitionTime)
      disp(['Time to run remaining loops ' num2str(SeqOut.SequenceTime*(SeqOut.Loops-Loop+1)+tb*(SeqOut.Loops-Loop+1),'% 10.1f') ' sec.']);
    else
      disp(['Time to run remaining loops ' num2str(SeqOut.LoopsRepetitionTime*(SeqOut.Loops-Loop+1)+SeqOut.SequenceTime*double(SeqOut.Loops==Loop),'% 10.1f') ' sec.']);
    end
  end
  clear tb
  % Run sequence
  SeqOut.PreProcessSequence = 0;
  SeqOut.StartSequence = 1;
  SeqOut.PollPPGfast = 1;
  SeqOut.GetRawData = 1;
  SeqOut.PostProcessSequence = 1;

  [~, SeqOut, data, data_1D] = set_sequence(HW, SeqOut, AQ, TX, Grad);
  SeqOut.data_1D = data_1D;
  % talker.mySequency.exctractArrayToFile(talker.mySequency.getCommandArray,'test.txt');

  if ~SeqOut.LoopsBreakExactly && ~isempty(SeqOut.LoopsBreak)
    init.Seq.StartSequenceTime = ...
      SeqOut.StartSequenceTime + SeqOut.SequenceTime + SeqOut.LoopsBreak;
  end

  if SeqOut.LoopsBreakExactly
    init.Seq.EndTimeFPGA = SeqOut.EndTimeFPGA + SeqOut.TR_Error/SeqOut.HW.MMRT(iDevice).fSystem*4;
    init.Seq.LoopCountStart = SeqOut.LoopCountEnd + SeqOut.TR_Error;
    init.Seq.TimeToNextSequence = SeqOut.TimeFromLastSequence;
    init.Seq.TimeFromLastSequence = SeqOut.TimeToNextSequence;
  end


  tr = 1;
  Seq.plotSeq = [];
  Seq.plotSeqTR = [];

  if Seq.LoopsRepetitionTime <= 30
    Seq.Reinitialize = 0;
    Seq.IgnoreTimingError = 1;
    Seq.TimeFromLastSequence = SeqOut.LoopsRepetitionTime - sum(SeqOut.SequenceTime);
    Seq.TimeToNextSequence = Seq.TimeFromLastSequence;
  else
    Seq.StartSequenceTime = SeqOut.StartSequenceTime + SeqOut.LoopsRepetitionTime;
    HW.tRepInit = 0.5;
  end

  iAQ = find([SeqOut.AQ(:).Channel] == Channel & [SeqOut.AQ(:).Device] == iDevice, 1, 'first');

  if Loop == 1
    SeqLoop = SeqOut;
    SeqLoop.data = data(iAQ);
    SeqLoop.loopdata.fft1_data = zeros([size(data(iAQ).fft1_data ,1), size(data(iAQ).fft1_data,2), size(data(iAQ).fft1_data,3), Seq.Loops]);
    SeqLoop.loopdata.data = zeros([size(data(iAQ).data,1), size(data(iAQ).data,2), size(data(iAQ).data,3), Seq.Loops]);
    SeqLoop.loopdata.time_all = zeros([size(data(iAQ).data,1), size(data(iAQ).data,2), size(data(iAQ).data,3), Seq.Loops]);
    SeqLoop.loopdata.phaseOffset = zeros(Seq.Loops, 1);
    SeqLoop.loopdata.fOffset = zeros(Seq.Loops, 1);
    SeqLoop.loopdata.StartTime = zeros(Seq.Loops, 1);
    init.Seq.StartSequenceTimeOffset = SeqOut.StartSequenceTime;
  end
  if Seq.plot
    fh6 = figure(6);
    ax6(1) = subplot(2,1,1, 'Parent', fh6);
    plot(ax6(1), ...
      repmat(data(iAQ).time_all(:,1,tr),1,3), ...
      [abs(data(iAQ).data(:,1,tr)), real(data(iAQ).data(:,1,tr)), imag(data(iAQ).data(:,1,tr))]*data(iAQ).Amplitude2Uin(1)*1e6);
    title(ax6(1), 'Acquired signal');
    ylabel(ax6(1), ['Amplitude in ' char(181) 'V']);
    xlabel(ax6(1), 'Time in s');
    ax6(2)=subplot(2,1,2, 'Parent', fh6);
    plot(ax6(2),((data(iAQ).f_fft1_data(:,1,tr))/SeqOut.AQ(iAQ).Frequency(tr)-1)*1e6,abs(data(iAQ).fft1_data(:,1,tr)).*data(iAQ).Amplitude2Uin(1)*1e6);
    xlabel(ax6(2), 'Frequency in ppm'); %offset ppm
    ylabel(ax6(2), ['Amplitude in ' char(181) 'V']);
    title(ax6(2), 'FFT of acquired signal');
    xlim(ax6(2), [-10 10]);
    set(ax6(2), 'XDir', 'reverse');
  end
  if Loop >= 1
    SeqLoop.loopdata.StartTime(Loop) = sum([SeqOut.StartSequenceTime -init.Seq.StartSequenceTimeOffset]);
    SeqLoop.loopdata.fft1_data(:,:,:,Loop) = data(iAQ).fft1_data;
    SeqLoop.loopdata.time_all(:,:,:,Loop) = data(iAQ).time_all;
    SeqLoop.loopdata.data(:,:,:,Loop) = data(iAQ).data;
    temp = squeeze(data(iAQ).data(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr));
    SeqLoop.loopdata.phaseOffset(Loop) = ...
      get_MeanPhaseWeighted(temp(find(data(iAQ).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.PhaseOffsetTimeStart,1,'first')   :   find(data(iAQ).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.PhaseOffsetTimeStop,1,'first')));
    SeqLoop.loopdata.fOffset(Loop) = ...
      get_MeanPhaseDiffWeighted(temp(find(data(iAQ).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStart,1,'first')   :   find(data(iAQ).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStop,1,'first')))*SeqOut.AQ(iAQ).fSample(1)/2/pi;
    HW.fLarmor = HW.fLarmor - SeqLoop.loopdata.fOffset(Loop);
  else
    temp = squeeze(data.data(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr));
    HW.fLarmor = HW.fLarmor - ...
      get_MeanPhaseDiffWeighted(temp(find(data(iAQ).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStart,1,'first')   :   find(data(iAQ).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStop,1,'first')))*SeqOut.AQ(iAQ).fSample(1)/2/pi;
  end
end


end
