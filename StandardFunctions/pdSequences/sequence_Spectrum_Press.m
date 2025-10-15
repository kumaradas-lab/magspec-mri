function [SeqLoop, mySave] = sequence_Spectrum_Press(HW, Seq, AQ, TX, Grad, mySave)
%% Acquire a signal using a PRESS sequence
%
%   [SeqLoop, mySave] = sequence_Spectrum_Press(HW, Seq, AQ, TX, Grad, mySave)
%
% A PRESS (Point RESolved Sprectroscopy) sequence is a multi spin echo sequence
% which acquires the second echo. The excitation and refocusing pulses are slice
% selective in orientations that are perpendicular to each other.
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
%       CorrectFrequencyLoop
%             Boolean value that selects whether the TX and AQ frequency of the
%             next loop should be corrected with the frequency offset that has
%             been determined from the measurement signal of the current loop.
%             (Default: true)
%
%       CorrectFrequencyAQChannel
%             Index for the AQ channel that is used for the frequency tracking
%             (and correction). By default, 1 corresponds to the 1H frequency
%             (HW.fLarmor) and 2 corresponds to the X frequency (HW.fLarmorX).
%             (Default: 1)
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
%       sizePhaseSpoil
%             5 element vector with spoiler size for the following spoilers:
%               1:  spoiler in slice(1) ("thickness") direction aligned with
%                   slice(1) rephase gradient and slice(3) rephase gradient.
%                   (Default: Inf)
%               2:  spoiler in slice(2) ("thicknessInversion") direction aligned
%                   with slice(2) dephase gradient and slice(2) rephase
%                   gradient. (Default: Inf)
%               3:  spoiler in slice(3) ("thicknessInversion2") direction
%                   aligned with slice(3) dephase gradient and slice(3) rephase
%                   gradient. (Default: Inf)
%               4:  spoiler in slice(2) ("thicknessInversion") direction aligned
%                   with slice(1) rephase gradient and slice(3) rephase
%                   gradient. (Default: sizePhaseSpoil(1) )
%               5:  spoiler in slice(3) ("thicknessInversion2") direction
%                   aligned with slice(1) rephase gradient and slice(3) rephase
%                   gradient. (Default: sizePhaseSpoil(1) )
%
%       phaseCycling
%             Boolean value to enable phase cycling of the excitation pulse.
%             This is helpful to reduce the signal of spins excited by the
%             refocusing (inversion) pulse. This effectively multiplies the
%             number of acquired PRESS echoes.
%
%       phaseCycleSteps
%             Scalar integer with the number of steps in the phase cycle
%             (default: 2). The steps are evenly distributed for a complete
%             cycle of 360 degrees. If Seq.AQSlice(1).phaseCycling is switched
%             off, the value of this field is set to 1.
%
%     dualNuclear
%           Boolean value to select whether the measurement should be run in
%           "dual frequency" mode. If false, the excitation pulses and
%           acquisitions are done at the frequency corresponding to HW.GammaDef.
%           If true, the excitation pulses and acquisition windows are mixed
%           from two frequencies corresponding to HW.GammaDef and HW.GammaX.
%           (Default: false)
%
%     plotAverage
%           Boolean value to select whether the (naive) average of all steps up
%           to that point in time are plotted after each average step. (Note:
%           Instead of this "naive" average, a frequency corrected average
%           should be used for spectrum analysis. See: get_Spectrum)
%           (Default: false)
%
%     Function_Prepare_Measurement
%           A function handle to a prepare function with the following
%           signature:
%                     [HW, Seq, AQ, TX, Grad] = @(HW, Seq, AQ, TX, Grad)
%           It is executed after the read, phase, and slice pulses ("logical
%           blocks") are created with get_ReadParameter, get_PhaseParameter,
%           and get_SliceParameter and before the corresponding pulses are
%           added to the final pulse program (AQ, TX, and Grad). It can be used
%           to modify the sequence before it is eventually executed.
%           Also see the programming notes below.
%           Alternatively, this can be a cell array of function handles with the
%           same signature as described before. In that case, the prepare
%           functions are executed in the order as they appear in the cell
%           array.
%           (Default: [])
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
% PROGRAMMING NOTES:
%
%   The pulse program is constructed using logical blocks that correspond to
%   slice, read, and phase encoders.
%   The following functions translate the logical blocks to the actual pulse
%   program:
%     For slice units:      get_SliceParameter
%     For read units:       get_ReadParameter
%     For phase units:      get_PhaseParameter
%   Please, see the documentation of these functions for further details.
%
%   When manipulating the pulse program with Seq.Function_Prepare_Measurement,
%   the settings for these units can be changed. Currently, the following units
%   are used:
%     Seq.Read(1)
%       Readout starting at the center of the PRESS echo
%     Seq.Slice(1)
%       90 degrees excitation pulse with slice selection
%     Seq.Slice(2)
%       180 degrees inversion pulse with slice selection perpendicular to
%       Seq.Slice(1)
%     Seq.Slice(3)
%       180 degrees inversion pulse with slice selection perpendicular to
%       Seq.Slice(1) and Seq.Slice(2)
%     Seq.Phase(1)
%       Spoiler in slice(1) ("thickness") direction aligned with slice(1)
%       rephase gradient and slice(3) rephase gradient.
%     Seq.Phase(2)
%       Spoiler in slice(2) ("thicknessInversion") direction aligned with
%       slice(2) dephase gradient and slice(2) rephase gradient.
%     Seq.Phase(3)
%       Spoiler in slice(3) ("thicknessInversion2") direction aligned with
%       slice(3) dephase gradient and slice(3) rephase gradient.
%     Seq.Phase(4)
%        Spoiler in slice(2) ("thicknessInversion") direction aligned with
%        slice(1) rephase gradient and slice(3) rephase gradient.
%     Seq.Phase(5)
%        Spoiler in slice(3) ("thicknessInversion2") direction aligned with
%        slice(1) rephase gradient and slice(3) rephase gradient.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2024 Pure Devices GmbH, Wuerzburg, Germany
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

if isemptyfield(Seq, 'Loops'), Seq.Loops = 1; end  % number of loops for averages
if isemptyfield(Seq, 'PreLoops')
  % number of pre-shots that are discarded in averages
  Seq.PreLoops = 1;
end
if isemptyfield(Seq, 'plot')
  % show plots with data of each single average step
  Seq.plot = 0;
end
if isemptyfield(Seq, 'plotAverage')
  % show plots with averaged data at each step
  Seq.plotAverage = Seq.plot && (Seq.Loops>1);
end
if isemptyfield(Seq, 'fSample'), Seq.fSample = HW.FindShim.fSample; end  % HW.RX.fSample = 125e6 Hz
if isemptyfield(Seq, 'fSampleFID'), Seq.fSampleFID = Seq.fSample; end  % HW.RX.fSample = 125e6 Hz
if isemptyfield(Seq, 'LoopsRepetitionTime'), Seq.LoopsRepetitionTime = HW.FindFrequencyPause; end
if isemptyfield(Seq, 'Find_Frequency_interval'), Seq.Find_Frequency_interval = 0; end  % Find Frequency interval

if isemptyfield(Seq, 'Spectro'), Seq.Spectro = []; end
if isemptyfield(Seq.Spectro, 'fOffsetTimeStart'), Seq.Spectro.fOffsetTimeStart = 0.03; end
if isemptyfield(Seq.Spectro, 'fOffsetTimeStop'), Seq.Spectro.fOffsetTimeStop = 0.09; end
if isemptyfield(Seq.Spectro, 'PhaseOffsetTimeStart'), Seq.Spectro.PhaseOffsetTimeStart = 0.03; end
if isemptyfield(Seq.Spectro, 'PhaseOffsetTimeStop'), Seq.Spectro.PhaseOffsetTimeStop = 0.04; end
if isemptyfield(Seq.Spectro, 'CorrectFrequencyLoop')
  Seq.Spectro.CorrectFrequencyLoop = true;
end
if isemptyfield(Seq.Spectro, 'CorrectFrequencyAQChannel')
  % use GammaDef signal for frequency drift by default
  Seq.Spectro.CorrectFrequencyAQChannel = 1;
end

if nargin < 3
  % matrix containing the starting times (center of the first Sample –
  % 0.5/AQ.fSample) of each acquisition window for each TR; up to 510 AQ windows
  % are possible per TR
  AQ.Start = [];
  % matrix containing the number of samples to be acquired for each AQ window
  % for each TR
  AQ.nSamples = [];
  % matrix containing the sampling frequency of each acquisition window for each
  % TR
  AQ.fSample = [];
  % matrix containing the mixing frequency of each acquisition windows for each
  % TR
  AQ.Frequency = [];
  % matrix containing the phase offset for each acquisition window
  AQ.Phase = [];
  % if set to 1 the TX and RX phase is matched at the beginning of the TR
  AQ.ResetPhases = [];
  % relative gain of the acquisition path (1/AQ.Gain = max input Amplitude).
  % Use HW.RX.GainDef for the best noise figure.
  AQ.Gain = [];
  % if there is a ‘1’ in the array, the settings for the window of the prior TR
  % are used. (reduces data traffic between the PC and the MRI device)
  AQ.Repeat = [];
end

if nargin < 4
  TX.Channel = HW.TX.ChannelDef;  % use transmitting channel 1 or 2
  % matrix containing the starting time of the RF pulses for each TR column by
  % column; up to 510 RF pulses are possible per TR (column)
  TX.Start = [];
  % matrix containing transmitting frequency for each pulse of each TR
  TX.Frequency = [];
  % matrix containing the duration of the RF pulses for each TR column by
  % column; up to 510 RF pulses are possible per TR
  TX.Duration = [];
  % matrix containing the amplitude for each pulse of each TR
  TX.Amplitude = [];
  % matrix containing the phase offset for each pulse of each TR
  TX.Phase = [];

  % array of times to define, when the blanking signal is applied prior to the
  % RF pulse; can be adjusted every TR
  TX.BlankOffset = [];
  % array of times to define how long blanking stays active after a RF pulse;
  % can be adjusted every TR
  TX.BlankPostset = [];
  % if there is a ‘1’ in the array, the pulses of the prior TR are used.
  % (reduces data traffic between the PC and the MRI device)
  TX.Repeat = [];
end

iDevice = 1;  % FIXME: Support multiple MMRT devices
Channel = 1;  % FIXME: Add support for multiple AQ channels?

if nargin < 5
  % time the values in Grad(t).Amp with the corresponding index are set
  [Grad(1:HW.Grad(iDevice).n).Time] = deal([]);
  % amplitude of the gradient at the time Grad(t).Time (linearly interpolated)
  [Grad(1:HW.Grad(iDevice).n).Amp] = deal([]);
  % additional shim; magnet shim is already considered in HW.MagnetShim.
  % Caution: Using high values over a long time will damage the gradient coils
  %          and amplifiers!
  [Grad(1:HW.Grad(iDevice).n).Shim] = deal([]);
  % if there is a one in the array the gradients of the prior TR are used.
  % (reduces data traffic between the PC and the MRI device)
  [Grad(1:HW.Grad(iDevice).n).Repeat] = deal([]);
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
if isemptyfield(Seq.AQSlice, 'inversion2FlipAngleComposite')
  Seq.AQSlice(1).inversion2FlipAngleComposite = Seq.AQSlice(1).inversionFlipAngleComposite;
end
if isemptyfield(Seq.AQSlice, 'inversion2FlipPhaseComposite')
  Seq.AQSlice(1).inversion2FlipPhaseComposite = Seq.AQSlice(1).inversionFlipAngleComposite;
end
if isemptyfield(Seq.AQSlice, 'excitationPhase'), Seq.AQSlice(1).excitationPhase = 0; end
if isemptyfield(Seq.AQSlice, 'inversionPhase'), Seq.AQSlice(1).inversionPhase = 90; end
if isemptyfield(Seq.AQSlice, 'phaseCycling')
  % acquire "PRESS echo" with phase of the excitation pulse cycled
  Seq.AQSlice(1).phaseCycling = 0;
end
if isemptyfield(Seq.AQSlice, 'phaseCycleSteps')
  % number of physe cycle steps with the phase of the excitation pulse
  % incremented by 360/phaseCycleSteps degrees.
  Seq.AQSlice(1).phaseCycleSteps = 2;
end
if ~Seq.AQSlice(1).phaseCycling
  Seq.AQSlice(1).phaseCycleSteps = 1;
end
if isemptyfield(Seq.AQSlice, 'excitationPhaseIncrement')
  Seq.AQSlice(1).excitationPhaseIncrement = 360/Seq.AQSlice(1).phaseCycleSteps;
end
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

if isemptyfield(Seq, 'dualNuclear')
  Seq.dualNuclear = false;
end

% FIXME: Add support for multiple devices
if Seq.dualNuclear
  if Seq.Spectro.CorrectFrequencyAQChannel > 2
    error('PD:sequence_Spectrum_Press:InvalidCorrectFrequencyAQChannel', ...
      'Seq.Spectro.CorrectFrequencyAQChannel must be set to either 1 (HW.fLarmor) or 2 (HW.fLarmorX).');
  end
else
  if Seq.Spectro.CorrectFrequencyAQChannel > 2
    error('PD:sequence_Spectrum_Press:InvalidCorrectFrequencyAQChannel', ...
      'Seq.Spectro.CorrectFrequencyAQChannel must be set to 1.');
  end
end

if numel(Seq.AQSlice.sizePhaseSpoil) < 4
  Seq.AQSlice.sizePhaseSpoil(4) = Seq.AQSlice.sizePhaseSpoil(1);
end
if numel(Seq.AQSlice.sizePhaseSpoil) < 5
  Seq.AQSlice.sizePhaseSpoil(5) = Seq.AQSlice.sizePhaseSpoil(1);
end

if ~isfield(Seq, 'Function_Prepare_Measurement')
  Seq.Function_Prepare_Measurement = [];
end


if isscalar(Seq.tRep)
  Seq.tRep = [Seq.LoopsRepetitionTime*ones(1, Seq.AQSlice(1).phaseCycleSteps-1), Seq.tRep];
end

Seq.P90tReps = 1:Seq.AQSlice(1).phaseCycleSteps;

%% initial values for all loops
init.AQ = AQ;
init.TX = TX;
init.Seq = Seq;
init.Grad = Grad;
LoopCount = 0;

%% loops
for loop = [PreLoops, 1:Seq.Loops]
  LoopCount = LoopCount+1;
  AQ = init.AQ;
  TX = init.TX;
  Grad = init.Grad;
  Seq = init.Seq;
  Seq.Loop = loop;
  Seq.LoopNameCount = LoopCount;


  if Seq.Loops+Seq.PreLoops>1
    fprintf('loops to run %d (%.1f s)\n', ...
      Seq.Loops-loop+1, (Seq.Loops-loop+1)*Seq.LoopsRepetitionTime);
  end


  % default parameters
  Seq = set_EmptyField(Seq, 'PreProcessSequence', 1);
  Seq = set_EmptyField(Seq, 'StartSequence', 1);
  Seq = set_EmptyField(Seq, 'PollPPGfast', 1);
  Seq = set_EmptyField(Seq, 'GetRawData', 1);
  Seq = set_EmptyField(Seq, 'PostProcessSequence', 1);


  %% create pulse sequence "elements"
  % Readout starting at the center of the PRESS echo
  Seq.Read(1).HzPerPixelMin = 1/Seq.tFID;
  Seq.Read(1).CenterOfReadout = 0;
  Seq.Read(1).UseCoordinate = Seq.AQSlice(1).ReadCoordinate;
  Seq.Read(1).GradTimeIntegralOffset = 0;
  Seq.Read(1).UseAtRepetitionTime = Seq.P90tReps;
  Seq.Read(1).Phase = Seq.AQSlice(1).excitationPhase;
  Seq.Read(1).PhaseIncrement = Seq.AQSlice(1).excitationPhaseIncrement;
  Seq.Read(1).distance = 0;
  Seq.Read(1).GradSign = 1;
  Seq.Read(1).useAQSlice = 1;
  Seq.Read(1).StartWithKSpaceCenter = Seq.AQSlice(1).StartWithKSpaceCenter;
  if Seq.dualNuclear
    Seq.Read(1).GammaX = HW.GammaX;
  end


  % 90 degrees slice excitation
  Seq.Slice(1).Pulse.Function = Seq.AQSlice(1).excitationPulse;
  Seq.Slice(1).Pulse.FlipAngleComposite = Seq.AQSlice(1).excitationFlipAngleComposite;
  Seq.Slice(1).Pulse.FlipPhaseComposite = Seq.AQSlice(1).excitationFlipPhaseComposite;
  Seq.Slice(1).Pulse.FlipAngle = Seq.AQSlice(1).excitationFlipAngle;
  Seq.Slice(1).Thickness = Seq.AQSlice(1).thickness;
  Seq.Slice(1).CenterOfPulse = -4*Seq.tEcho/2;
  Seq.Slice(1).UseCoordinate = 1;
  Seq.Slice(1).UseAtRepetitionTime = Seq.Read(1).UseAtRepetitionTime;
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
  if Seq.dualNuclear
    Seq.Slice(1).GammaX = HW.GammaX;
  end


  % 180 degrees (slice) inversion
  Seq.Slice(2).Pulse.Function = Seq.AQSlice(1).inversionPulse2;
  Seq.Slice(2).Pulse.FlipAngleComposite = Seq.AQSlice(1).inversionFlipAngleComposite;
  Seq.Slice(2).Pulse.FlipPhaseComposite = Seq.AQSlice(1).inversionFlipPhaseComposite;
  Seq.Slice(2).Pulse.FlipAngle = Seq.AQSlice(1).inversionFlipAngle;
  Seq.Slice(2).Thickness = Seq.AQSlice(1).thicknessInversion;
  Seq.Slice(2).CenterOfPulse = -3*Seq.tEcho/2;
  Seq.Slice(2).UseCoordinate = 2;
  Seq.Slice(2).UseAtRepetitionTime = Seq.Read(1).UseAtRepetitionTime;
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
  if Seq.dualNuclear
    Seq.Slice(2).GammaX = HW.GammaX;
  end

  % 180 degrees (slice) inversion2
  Seq.Slice(3).Pulse.Function = Seq.AQSlice(1).inversionPulse;
  Seq.Slice(3).Pulse.FlipAngleComposite = Seq.AQSlice(1).inversionFlipAngleComposite;
  Seq.Slice(3).Pulse.FlipPhaseComposite = Seq.AQSlice(1).inversionFlipPhaseComposite;
  Seq.Slice(3).Pulse.FlipAngle = Seq.AQSlice(1).inversionFlipAngle;
  Seq.Slice(3).Thickness = Seq.AQSlice(1).thicknessInversion2;
  Seq.Slice(3).CenterOfPulse = -1*Seq.tEcho/2;
  Seq.Slice(3).UseCoordinate = 3;
  Seq.Slice(3).UseAtRepetitionTime = Seq.Read(1).UseAtRepetitionTime;
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
  if Seq.dualNuclear
    Seq.Slice(3).GammaX = HW.GammaX;
  end

  % generate readout timings
  Seq = get_SliceParameter(Seq, HW);

  % Seq.Slice(1).CenterOfRephase=Seq.Slice(2).CenterOfDephase;
  Seq.Slice(2).tEC = ...
    max(Seq.AQSlice(1).tCrusherSliceEqual/2-Seq.Slice(2).GradLength/2+Seq.Slice(2).tRamp, ...
        Seq.AQSlice(1).tCrusherEC);
  Seq.Slice(3).tEC = ...
    max(Seq.AQSlice(1).tCrusherSliceEqual/2-Seq.Slice(3).GradLength/2+Seq.Slice(3).tRamp, ...
        Seq.AQSlice(1).tCrusherEC);
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


  % spoiler in slice "thickness" direction aligned with
  % Seq.Slice(1).GradRephase and Seq.Slice(3).GradRephase
  Seq.Phase(1).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(1);
  Seq.Phase(1).nPhase = 1;
  Seq.Phase(1).PhaseOS = 2;
  Seq.Phase(1).StepOrder = [1, 1];
  Seq.Phase(1).CenterOfDephase = Seq.Slice(1).CenterOfRephase;
  Seq.Phase(1).CenterOfRephase = Seq.Slice(3).CenterOfRephase;
  Seq.Phase(1).GradDephaseLength = Seq.Slice(1).GradRephaseLength;
  Seq.Phase(1).GradRephaseLength = Seq.Slice(3).GradRephaseLength;
  Seq.Phase(1).UseCoordinate = 1;
  Seq.Phase(1).UseAtRepetitionTime = Seq.Read(1).UseAtRepetitionTime;
  Seq.Phase(1).GradDephaseSign = -1;
  Seq.Phase(1).GradRephaseSign = 1;

  % spoiler in slice "thicknessInversion" direction aligned with
  % Seq.Slice(2).CenterOfDephase and Seq.Slice(2).CenterOfRephase
  Seq.Phase(2).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(2);
  Seq.Phase(2).nPhase = 1;
  Seq.Phase(2).PhaseOS = 2;
  Seq.Phase(2).StepOrder = [1, 1];
  Seq.Phase(2).CenterOfDephase = Seq.Slice(2).CenterOfDephase;
  Seq.Phase(2).CenterOfRephase = Seq.Slice(2).CenterOfRephase;
  Seq.Phase(2).GradDephaseLength = Seq.Slice(2).GradDephaseLength;
  Seq.Phase(2).GradRephaseLength = Seq.Slice(2).GradRephaseLength;
  Seq.Phase(2).UseCoordinate = 2;
  Seq.Phase(2).UseAtRepetitionTime = Seq.Read(1).UseAtRepetitionTime;
  Seq.Phase(2).GradDephaseSign = -1;
  Seq.Phase(2).GradRephaseSign = -1;

  % spoiler in slice "thicknessInversion2" direction aligned with
  % Seq.Slice(3).CenterOfDephase and Seq.Slice(3).CenterOfRephase
  Seq.Phase(3).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(3);
  Seq.Phase(3).nPhase = 1;
  Seq.Phase(3).PhaseOS = 2;
  Seq.Phase(3).StepOrder = [1, 1];
  Seq.Phase(3).CenterOfDephase = Seq.Slice(3).CenterOfDephase;
  Seq.Phase(3).CenterOfRephase = Seq.Slice(3).CenterOfRephase;
  Seq.Phase(3).GradDephaseLength = Seq.Slice(3).GradDephaseLength;
  Seq.Phase(3).GradRephaseLength = Seq.Slice(3).GradRephaseLength;
  Seq.Phase(3).UseCoordinate = 3;
  Seq.Phase(3).UseAtRepetitionTime = Seq.Read(1).UseAtRepetitionTime;
  Seq.Phase(3).GradDephaseSign = -1;
  Seq.Phase(3).GradRephaseSign = -1;

  % spoiler in slice "thicknessInversion" direction aligned with
  % Seq.Slice(1).GradRephase and Seq.Slice(3).GradRephase
  Seq.Phase(4).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(4);
  Seq.Phase(4).nPhase = 1;
  Seq.Phase(4).PhaseOS = 2;
  Seq.Phase(4).StepOrder = [1, 1];
  Seq.Phase(4).CenterOfDephase = Seq.Slice(1).CenterOfRephase;
  Seq.Phase(4).CenterOfRephase = Seq.Slice(3).CenterOfRephase;
  Seq.Phase(4).GradDephaseLength = Seq.Slice(1).GradRephaseLength;
  Seq.Phase(4).GradRephaseLength = Seq.Slice(3).GradRephaseLength;
  Seq.Phase(4).UseCoordinate = 2;
  Seq.Phase(4).UseAtRepetitionTime = Seq.Read(1).UseAtRepetitionTime;
  Seq.Phase(4).GradDephaseSign = 1;
  Seq.Phase(4).GradRephaseSign = -1;

  % spoiler in slice "thicknessInversion2" direction aligned with
  % Seq.Slice(1).GradRephase and Seq.Slice(3).GradRephase
  Seq.Phase(5).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(5);
  Seq.Phase(5).nPhase = 1;
  Seq.Phase(5).PhaseOS = 2;
  Seq.Phase(5).StepOrder = [1, 1];
  Seq.Phase(5).CenterOfDephase = Seq.Slice(1).CenterOfRephase;
  Seq.Phase(5).CenterOfRephase = Seq.Slice(3).CenterOfRephase;
  Seq.Phase(5).GradDephaseLength = Seq.Slice(1).GradRephaseLength;
  Seq.Phase(5).GradRephaseLength = Seq.Slice(3).GradRephaseLength;
  Seq.Phase(5).UseCoordinate = 3;
  Seq.Phase(5).UseAtRepetitionTime = Seq.Read(1).UseAtRepetitionTime;
  Seq.Phase(5).GradDephaseSign = -1;
  Seq.Phase(5).GradRephaseSign = 1;


  % generate readout timings
  Seq = get_ReadParameter(Seq, HW);
  % generate slice parameters
  Seq = get_SliceParameter(Seq, HW);
  % generate phase parameters
  Seq = get_PhaseParameter(Seq, HW);
  % Gradients of image encoding


  %% User function hook for sequence manipulation
  if ~isempty(Seq.Function_Prepare_Measurement)
    if iscell(Seq.Function_Prepare_Measurement)
      for iPrep = 1:numel(Seq.Function_Prepare_Measurement)
        [HW, Seq, AQ, TX, Grad] = Seq.Function_Prepare_Measurement{iPrep}(HW, Seq, AQ, TX, Grad);
      end
    else
      [HW, Seq, AQ, TX, Grad] = Seq.Function_Prepare_Measurement(HW, Seq, AQ, TX, Grad);
    end
  end


  %% merge all elements into the complete pulse program
  % AQ
  AQ = add_AQ(AQ, Seq.Read(1).AQ);

  % TX
  for t = 1:3
    TX = add_TX(TX, Seq.Slice(t).TX);
  end

  % Grad
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

  if Seq.AQSlice(1).sizePhaseSpoil(4) < 1
    Grad = add_Grad(Grad, Seq.Phase(4).GradDephase);
    Grad = add_Grad(Grad, Seq.Phase(4).GradRephase);
  end

  if Seq.AQSlice(1).sizePhaseSpoil(5) < 1
    Grad = add_Grad(Grad, Seq.Phase(5).GradDephase);
    Grad = add_Grad(Grad, Seq.Phase(5).GradRephase);
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
      % PreProcessSequence sequence
      [~, SeqOut] = set_sequence(HW, Seq, AQ, TX, Grad);
    else
      % PreProcessSequence sequence only
      [~, SeqLoop] = set_sequence(HW, Seq, AQ, TX, Grad);
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
    elseif (SeqOut.LoopNameCount >= 2) ...
        && (SeqOut.LoopNameCount ~= numel(SeqOut.LoopName))
      SeqOut.Reinitialize = 0;
      SeqOut.TimeFromLastSequence = SeqOut.LoopsBreak+SeqOut.tOffset(1);
      SeqOut.TimeToNextSequence = SeqOut.LoopsBreak+SeqOut.tOffset(1);
    elseif SeqOut.LoopNameCount == numel(SeqOut.LoopName)
      SeqOut.Reinitialize = 0;
      SeqOut.TimeFromLastSequence = SeqOut.LoopsBreak+SeqOut.tOffset(1);
      SeqOut.TimeToNextSequence = [];
    end
  end

  if isempty(SeqOut.LoopsBreak)
    tb = 1;
  else
    tb = SeqOut.LoopsBreak;
  end
  if strcmp(SeqOut.LoopName{SeqOut.LoopNameCount}, 'normal') && (SeqOut.Loops > 1)
    if isempty(SeqOut.LoopsRepetitionTime)
      fprintf('Time to run remaining loops % 10.1f sec.\n', ...
        SeqOut.SequenceTime*(SeqOut.Loops-loop+1)+tb*(SeqOut.Loops-loop+1));
    else
      fprintf('Time to run remaining loops % 10.1f sec.\n', ...
        SeqOut.LoopsRepetitionTime*(SeqOut.Loops-loop+1)+SeqOut.SequenceTime*double(SeqOut.Loops==loop));
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
    init.Seq.EndTimeFPGA = SeqOut.EndTimeFPGA ...
      + SeqOut.TR_Error ./ [SeqOut.HW.MMRT(:).fSystem] * 4;
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
    Seq.TimeFromLastSequence = SeqOut.LoopsRepetitionTime*numel(Seq.tRep) - sum(SeqOut.SequenceTime);
    Seq.TimeToNextSequence = Seq.TimeFromLastSequence;
  else
    Seq.StartSequenceTime = SeqOut.StartSequenceTime + SeqOut.LoopsRepetitionTime*numel(Seq.tRep);
    HW.tRepInit = 0.5;
  end

  iAQ = find([SeqOut.AQ(:).Channel] == Channel & [SeqOut.AQ(:).Device] == iDevice, 1, 'first');
  iData = find([data(:).channel] == Channel & [data(:).device] == iDevice);

  if loop == 1
    SeqLoop = SeqOut;
    SeqLoop.data = data(iAQ);
    SeqLoop.loopdata.fft1_data = ...
      zeros([size(data(iAQ).fft1_data ,1), size(data(iAQ).fft1_data,2), size(data(iAQ).fft1_data,3), Seq.Loops]);
    SeqLoop.loopdata.data = ...
      zeros([size(data(iAQ).data,1), size(data(iAQ).data,2), size(data(iAQ).data,3), Seq.Loops]);
    SeqLoop.loopdata.time_all = ...
      zeros([size(data(iAQ).data,1), size(data(iAQ).data,2), size(data(iAQ).data,3), Seq.Loops]);
    for iAQX = 1:numel(iData)
      SeqLoop.loopdata.f_fft1_data(:,1,:,1,iAQX) = data(iData(iAQX)).f_fft1_data(:,1,:);
      SeqLoop.loopdata.Amplitude2Uin(iAQX) = data(iData(iAQX)).Amplitude2Uin(1);
    end
    SeqLoop.loopdata.phaseOffset = zeros(Seq.Loops, 1);
    SeqLoop.loopdata.fOffset = zeros(Seq.Loops, 1);
    SeqLoop.loopdata.StartTime = zeros(Seq.Loops, 1);
    init.Seq.StartSequenceTimeOffset = SeqOut.StartSequenceTime;
  end
  if Seq.plot
    for iAQX = 1:numel(iData)
      h_fig_data = figure(5+iAQX);
      h_ax_data(1) = subplot(2,1,1, 'Parent', h_fig_data);
      % for these "intermediate" plot just take the average of the phase cycling
      % steps (without phase correction)
      plot(h_ax_data(1), ...
        repmat(data(iData(iAQX)).time_all(:,1,tr),1,3), ...
        [abs(mean(data(iData(iAQX)).data(:,1,:), 3)), ...
         real(mean(data(iData(iAQX)).data(:,1,:), 3)), ...
         imag(mean(data(iData(iAQX)).data(:,1,:), 3))] * data(iData(iAQX)).Amplitude2Uin(1) * 1e6);
      title(h_ax_data(1), 'Acquired signal');
      ylabel(h_ax_data(1), ['Amplitude in ' char(181) 'V']);
      xlabel(h_ax_data(1), 'Time in s');
      grid(h_ax_data(1), 'on');
      h_ax_data(2)=subplot(2,1,2, 'Parent', h_fig_data);
      if iAQX == iAQ
        freq = SeqOut.AQ(iAQ).Frequency(tr);
      else
        freq = SeqOut.AQ(iAQ).FrequencyX(tr);
      end
      plot(h_ax_data(2), ...
        ((data(iData(iAQX)).f_fft1_data(:,1,tr))/freq-1) * 1e6, ...
        abs(mean(data(iData(iAQX)).fft1_data(:,1,:), 3)) ...
        .* data(iData(iAQX)).Amplitude2Uin(1) * 1e6);
      xlabel(h_ax_data(2), 'Chemical shift in ppm');  % offset ppm
      ylabel(h_ax_data(2), ['Amplitude in ' char(181) 'V']);
      title(h_ax_data(2), 'FFT of acquired signal');
      xlim(h_ax_data(2), [-10, 10]);
      set(h_ax_data(2), 'XDir', 'reverse');
      grid(h_ax_data(2), 'on');
      drawnow();
    end
  end

  iAQ_ref = iData(Seq.Spectro.CorrectFrequencyAQChannel);  % AQ mixer/channel for frequency tracking
  if loop >= 1
    SeqLoop.loopdata.StartTime(loop) = sum([SeqOut.StartSequenceTime, -init.Seq.StartSequenceTimeOffset]);
    for iAQX = 1:numel(iData)
      SeqLoop.loopdata.fft1_data(:,:,:,loop,iAQX) = data(iData(iAQX)).fft1_data;
      SeqLoop.loopdata.time_all(:,:,:,loop,iAQX) = data(iData(iAQX)).time_all;
      SeqLoop.loopdata.data(:,:,:,loop,iAQX) = data(iData(iAQX)).data;
    end
    temp = squeeze(data(iData(iAQ_ref)).data(1:SeqOut.AQ(iAQ).nSamples(tr),1,end));
    SeqLoop.loopdata.phaseOffset(loop) = ...
      get_MeanPhaseWeighted(...
        temp( find(data(iData(iAQ_ref)).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.PhaseOffsetTimeStart,1,'first') ...
             :find(data(iData(iAQ_ref)).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.PhaseOffsetTimeStop,1,'first')));
    SeqLoop.loopdata.fOffset(loop) = ...
      get_MeanPhaseDiffWeighted(...
        temp( find(data(iData(iAQ_ref)).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStart,1,'first') ...
             :find(data(iData(iAQ_ref)).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStop,1,'first'))) ...
      * SeqOut.AQ(iAQ).fSample(1)/2/pi;
    if Seq.Spectro.CorrectFrequencyAQChannel == 2
      % "scale" frequency offset at X frequency to 1H frequency
      SeqLoop.loopdata.fOffset(loop) = SeqLoop.loopdata.fOffset(loop) ...
        / HW.GammaX * HW.GammaDef;
      % FIXME: Should we do something with SeqLoop.loopdata.phaseOffset?
    end
    if Seq.Spectro.CorrectFrequencyLoop
      HW.fLarmor = HW.fLarmor - SeqLoop.loopdata.fOffset(loop);
    end

    if Seq.plotAverage
      % plot "naive" average with the measurements up to this point
      for iAQX = 1:numel(iData)
        h_fig_avg_data = figure(15+iAQX); clf(h_fig_avg_data);
        h_ax_avg_data(1) = subplot(2,1,1, 'Parent', h_fig_avg_data);
        % naive avaraging (without independent frequency shifts)
        avgData = mean(mean(SeqLoop.loopdata.data(:,:,:,1:loop,iAQX), 3), 4);
        plot(h_ax_avg_data(1),...
          repmat(data(iData(iAQX)).time_all(:,1,tr), 1, 3), ...
          [abs(avgData(:,1,tr)), real(avgData(:,1,tr)), imag(avgData(:,1,tr))]*data(iAQ).Amplitude2Uin(1)*1e6);
        title(h_ax_avg_data(1), 'Averaged signal');
        ylabel(h_ax_avg_data(1), ['Amplitude in ' char(181) 'V']);
        xlabel(h_ax_avg_data(1), 'Time in s');
        grid(h_ax_avg_data(1), 'on');
        h_ax_avg_data(2) = subplot(2,1,2, 'Parent', h_fig_avg_data);
        avgData_fft = fftshift(ifft(avgData(:,1,tr), [], 1));
        if iAQX == 1
          freq = SeqOut.AQ(iAQ).Frequency(tr);
        else
          freq = SeqOut.AQ(iAQ).FrequencyX(tr);
        end
        plot(h_ax_avg_data(2), ...
          ((data(iData(iAQX)).f_fft1_data(:,1,tr))/freq-1)*1e6, ...
          abs(avgData_fft).*data(iData(iAQX)).Amplitude2Uin(1)*1e6);
        xlabel(h_ax_avg_data(2), 'Frequency in ppm');  % offset ppm
        ylabel(h_ax_avg_data(2), ['Amplitude in ' char(181) 'V']);
        title(h_ax_avg_data(2), 'FFT of averaged signal');
        xlim(h_ax_avg_data(2), [-10, 10]);
        set(h_ax_avg_data(2), 'XDir', 'reverse');
        grid(h_ax_avg_data(2), 'on');
        drawnow();
      end
    end

  elseif Seq.Spectro.CorrectFrequencyLoop
    temp = squeeze(data(iData(iAQ_ref)).data(1:SeqOut.AQ(iAQ).nSamples(tr),1,end));
    fOffset = ...
      get_MeanPhaseDiffWeighted(...
        temp( find(data(iData(iAQ_ref)).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStart,1,'first') ...
             :find(data(iData(iAQ_ref)).time_all(1:SeqOut.AQ(iAQ).nSamples(tr),1,tr)>SeqOut.Spectro.fOffsetTimeStop,1,'first'))) ...
      * SeqOut.AQ(iAQ).fSample(1)/2/pi;
    HW.fLarmor = HW.fLarmor - double(fOffset);
  end
end

end
