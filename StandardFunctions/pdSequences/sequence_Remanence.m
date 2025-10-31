function SeqOut = sequence_Remanence(HW, Seq, AQ, TX, Grad)
%% "Imprint" remanence in pole shoes and probe
%
%   SeqOut = sequence_Remanence(HW, Seq, AQ, TX, Grad)
%
% A series of alternating gradient pulses with decreasing amplitude is emitted
% on the selected gradient coils to "imprint" a more or less reproducible
% remanence pattern in pole shoes and probe.
%
% This sequence can be run before other sequences that depend on remanence.
%
% INPUT:
%
%   HW
%       HW object or structure
%
%   Seq
%       Structure with the following (optional) field. If the fields are omitted
%       or empty, default values are used.
%
%     RemanenceGradients
%         Logic indexing vector with the gradients that are included in the
%         remanence pulse program. (Default: HW.Grad.ShimGradients)
%
%     GradLength
%         The duration of each gradient pulse in seconds. (Default: 0.5e-3)
%
%     tRamp
%         The ramp time for the gradient pulses in seconds. (Default:
%         HW.Grad.tRamp)
%
%     timeStep
%         The time for each gradient pulse in seconds. (Default:
%         Seq.GradLength - Seq.tRamp)
%
%     nSteps
%         The number of steps of decreasing amplitude. (Default: 32)
%
%     tRep
%         Repetition time for the entire pulse program in seconds. (Default:
%         sum(Seq.RemanenceGradients)*2*Seq.timeStep*Seq.nSteps)
%
%     maxAmpRelative
%         A factor to the maximum gradient amplitude (HW.Grad.MaxAmp) at the
%         respective channel. This can be a scalar or a vector. If it is a
%         vector, it must match the size of Seq.RemanenceGradients. (Default:
%         0.9)
%
%     iDevice
%         Index for the used MMRT device. (Default: 1)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% default values
if nargin < 2, Seq = struct(); end
if nargin < 3, AQ = struct(); end
if nargin < 4, TX = struct(); end
if nargin < 5, Grad = struct(); end

if isemptyfield(Seq, 'iDevice'), Seq.iDevice = 1; end

if isemptyfield(Seq, 'RemanenceGradients')
  Seq.RemanenceGradients = HW.Grad(Seq.iDevice).ShimGradients;
end

if isemptyfield(Seq, 'GradLength'), Seq.GradLength = 0.5e-3; end
if isemptyfield(Seq, 'tRamp'), Seq.tRamp = HW.Grad(Seq.iDevice).tRamp; end
if isemptyfield(Seq, 'timeStep')
  % Seq.timeStep = Seq.GradLength + HW.Grad(Seq.iDevice).tEC;
  Seq.timeStep = Seq.GradLength - Seq.tRamp;
end
if isemptyfield(Seq, 'nSteps'), Seq.nSteps = 32; end

% minimum tRep
if isemptyfield(Seq, 'tRep')
  Seq.tRep = sum(Seq.RemanenceGradients)*2*Seq.timeStep*Seq.nSteps;
end

if isemptyfield(Seq, 'maxAmpRelative'), Seq.maxAmpRelative = 0.9; end


if isemptyfield(Grad(1), 'Time'), Grad(1).Time = []; end
if isemptyfield(Grad(1), 'Amp'), Grad(1).Amp = []; end


%% check input
if (Seq.GradLength/2 - Seq.tRamp/2 > Seq.timeStep) || ...
    ((sum(Seq.RemanenceGradients)==1) && (Seq.GradLength/2 > Seq.timeStep))
  % FIXME: Do we need to enforce this?
  error('PD:sequence_Remanence:TooFast', ...
    'Seq.timeStep is too fast.');
end

if isscalar(Seq.maxAmpRelative)
  Seq.maxAmpRelative = Seq.maxAmpRelative * ones(size(Seq.RemanenceGradients));
end

if size(Seq.maxAmpRelative) ~= size(Seq.RemanenceGradients)
  error('PD:sequence_Remanence:WrongSizeMaxAmpRelative', ...
    'Seq.maxAmpRelative must have a compatible size.');
end

%% prepare pulse program
useGrad = find(Seq.RemanenceGradients);

for iRem = 1:numel(useGrad)
  Seq.Phase(iRem).nPhase = Seq.nSteps;
  Seq.Phase(iRem).PhaseOS = 2;  % only take steps until "k-space center"
  nSteps = Seq.Phase(iRem).nPhase*Seq.Phase(iRem).PhaseOS;
  % Each block must be in a separate tRep for the generating function.
  Seq.Phase(iRem).UseAtRepetitionTime = 1:nSteps;

  Seq.Phase(iRem).GradDephaseLength = Seq.GradLength;
  Seq.Phase(iRem).GradRephaseLength = Seq.GradLength;

  Seq.Phase(iRem).CenterOfDephase = (0:nSteps-1) * 2*numel(useGrad)*Seq.timeStep + ...
    (iRem-1)*Seq.timeStep;
  Seq.Phase(iRem).CenterOfRephase = Seq.Phase(iRem).CenterOfDephase + ...
    numel(useGrad)*Seq.timeStep;

  maxGradTimeIntegral = Seq.maxAmpRelative(iRem) * ...
    HW.Grad(Seq.iDevice).MaxAmp(useGrad(iRem)) * ...
    (Seq.GradLength - Seq.tRamp);

  Seq.Phase(iRem).Gamma = HW.GammaDef;

  % calculate size to have selected gradient amplitude
  Seq.Phase(iRem).sizePhase = Seq.Phase(1).nPhase * pi / Seq.Phase(1).Gamma / ...
    maxGradTimeIntegral;
  Seq.Phase(iRem).UseCoordinate = useGrad(iRem);

  Seq.Phase(iRem).alfa = 0;
  Seq.Phase(iRem).phi = 0;
  Seq.Phase(iRem).theta = 0;
end

Seq = get_PhaseParameter(Seq, HW);

% fold everything (until "k-space center") into a single tRep
for iRem = 1:numel(useGrad)
  for iGr = numel(Seq.Phase(1).GradDephase):-1:1
    GradDephase(iGr).Amp = reshape(Seq.Phase(iRem).GradDephase(iGr).Amp(:,1:Seq.nSteps), [], 1);
    GradDephase(iGr).Time = reshape(Seq.Phase(iRem).GradDephase(iGr).Time(:,1:Seq.nSteps), [], 1);
    GradRephase(iGr).Amp = reshape(Seq.Phase(iRem).GradRephase(iGr).Amp(:,1:Seq.nSteps), [], 1);
    GradRephase(iGr).Time = reshape(Seq.Phase(iRem).GradRephase(iGr).Time(:,1:Seq.nSteps), [], 1);
  end

  Grad = add_Grad(Grad, GradDephase);
  Grad = add_Grad(Grad, GradRephase);
end


%% run pulse program
[~, SeqOut] = set_sequence(HW, Seq, AQ, TX, Grad);


end
