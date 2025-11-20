function [HW, Seq, AQ, TX, Grad] = Prepare_Flash_EC(HW, Seq, AQ, TX, Grad)
%% Prepare function that adds gradient pulses to measure eddy current effects
%
%   [HW, Seq, AQ, TX, Grad] = Prepare_Flash_EC(HW, Seq, AQ, TX, Grad)
%
%
% INPUT:
%
%   The Seq structure may contain the field "EC" which is a structure with the
%   following optional fields:
%
%     ampGrad
%       Amplitude of the eddy current gradient pulse in Tesla. (Default: 10e-3)
%
%     useCoordinate
%       Coordinate for the eddy current gradient pulse (1, 2, or 3).
%       (Default: 1)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2022-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% check input
if ~isfield(Seq, 'EC'),  Seq.EC = struct();  end

if isemptyfield(Seq.EC, 'ampGrad'),  Seq.EC.ampGrad = 10e-3;  end
if isemptyfield(Seq.EC, 'useCoordinate'),  Seq.EC.useCoordinate = 1;  end


%% adjust pulse program
% end eddy current pulse with start of acquisition window
tDephaseOffset = Seq.Read(1).GradDephaseLength - Seq.Read(1).tEC;

GradTimeIntegral = (Seq.Read(1).GradDephaseLength - Seq.Read(1).tRamp) * Seq.EC.ampGrad;

% add gradients that cause an eddy current immediately before the read out
nPhase = numel(Seq.Phase);
iPhase = nPhase+1;
Seq.EC.iPhase = iPhase;  % informationally
Seq.Phase(iPhase).sizePhase = pi/Seq.Phase(1).Gamma / GradTimeIntegral;
Seq.Phase(iPhase).nPhase = 1;
Seq.Phase(iPhase).PhaseOS = 2;
Seq.Phase(iPhase).StepOrder = [1, 1];
% Seq.Phase(iPhase).GradDephaseSign = ones(1, numel(Seq.Phase(4).UseAtRepetitionTime));
% Seq.Phase(iPhase).GradDephaseSign(2:2:end) = -1;
% Seq.Phase(iPhase).GradRephaseSign = -Seq.Phase(iPhase).GradDephaseSign;
Seq.Phase(iPhase).CenterOfDephase = Seq.Read(1).CenterOfDephase - tDephaseOffset;
Seq.Phase(iPhase).CenterOfRephase = Seq.Read(1).CenterOfDephase + Seq.Read(1).tEC;
Seq.Phase(iPhase).GradDephaseLength = Seq.Read(1).GradDephaseLength;
Seq.Phase(iPhase).GradRephaseLength = Seq.Read(1).GradDephaseLength;
Seq.Phase(iPhase).UseCoordinate = Seq.EC.useCoordinate;
Seq.Phase(iPhase).GradTimeDelay = Seq.Phase(4).GradTimeDelay;
% Seq.Phase(iPhase).UseAtRepetitionTime = Seq.Phase(4).UseAtRepetitionTime;
Seq.Phase(iPhase).UseAtRepetitionTime = Seq.Phase(1).UseAtRepetitionTime;
% Seq.Phase(iPhase).tRamp = HW.Grad(Seq.AQSlice(1).iDevice).tRamp;  % FIXME: Is it possible to only change tRamp of the rephase pulse?


% move encoders away from the read out to make room for the eddy current pulse
for iPhase = 1:min(6, numel(Seq.Phase))
  Seq.Phase(iPhase).CenterOfDephase = Seq.Phase(iPhase).CenterOfDephase - tDephaseOffset;
end

for iRead = 1:min(2, numel(Seq.Read))
  Seq.Read(iRead).CenterOfDephase = Seq.Read(iRead).CenterOfDephase - tDephaseOffset;
end


%% (re-)create pulse program for adjusted units
Seq = get_PhaseParameter(Seq, HW);
Seq = get_ReadParameter(Seq, HW);


%% add gradients to pulse program
Grad = add_Grad(Grad, Seq.Phase(end).GradDephase);
Grad = add_Grad(Grad, Seq.Phase(end).GradRephase);


end
