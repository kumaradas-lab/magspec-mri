function [HW, Seq, AQ, TX, Grad] = prepare_T2T2(HW, Seq, AQ, TX, Grad)
%% Prepare CPMG sequence using a CPMG echo train
%
%     [HW, Seq, AQ, TX, Grad] = prepare_T2T2(HW, Seq, AQ, TX, Grad)
%
% Add CPMG echo trains as a preparation block for T2-T2 map acquisition.
%
% Important input parameters:
%
%   Seq.T2Prepare
%     Structure with the following (optional) fields:
%
%       tEcho
%         Echo time in the preparation CPMG echo train in seconds.
%         (Default: Seq.tEcho)
%
%       tEchoTrain
%         (Approximate) total time of the preparation CPMG echo train in
%         seconds. nEcho takes precedence.
%         (Default: Seq.tEchoTrain)
%
%       nEcho
%         Number of echoes in preparation CPMG echo train in seconds.
%         (Default: round(Seq.T2Prepare.tEchoTrain/Seq.T2Prepare.tEcho))
%
%       excitationPulse
%         Handle to a pulse shape function for the 90 degrees excitation pulse
%         for the preparation CPMG echo train (default: Seq.excitationPulse).
%
%       refocusingPulse
%         Handle to a pulse shape function for the 180 degrees refocusing pulse
%         for the preparation CPMG echo train (default: Seq.refocusingPulse).
%
%
% See also: sequence_RecoveryCPMG
%
% ------------------------------------------------------------------------------
% (C) Copyright 2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% Default parameters
iDevice = 1;  % FIXME: Add support for multiple devices

if isemptyfield(Seq, 'T2Prepare'),  Seq.T2Prepare = struct();  end

if isemptyfield(Seq.T2Prepare, 'tEcho'), Seq.T2Prepare.tEcho = Seq.tEcho; end
if isemptyfield(Seq.T2Prepare, 'tEchoTrain'), Seq.T2Prepare.tEchoTrain = Seq.tEchoTrain; end
if isemptyfield(Seq.T2Prepare, 'nEcho')
  Seq.T2Prepare.nEcho = round(Seq.T2Prepare.tEchoTrain/Seq.T2Prepare.tEcho);
end
if isemptyfield(Seq.T2Prepare, 'excitationPulse')
  Seq.T2Prepare.excitationPulse = Seq.excitationPulse;
end
if isemptyfield(Seq.T2Prepare, 'refocusingPulse')
  Seq.T2Prepare.refocusingPulse = Seq.refocusingPulse;
end
if isemptyfield(Seq.T2Prepare, 'Phase')
  Seq.T2Prepare.Phase = 180;
end


%% Add preparation CPMG echo train

% elements are actions x tReps x nTau1
nTau1 = size(TX.Start, 3);  % actual number of tau1 and steady state repetitions

%% prepend tReps for preparation
Seq.T2Prepare.tRepsPrepare = Seq.T2Prepare.nEcho+2;
Seq.tRep = cat(2, ...
  repmat([0.5, ones(1,Seq.T2Prepare.nEcho-1), 0.5]*Seq.T2Prepare.tEcho, [1,1,nTau1]), ...
  Seq.tRep);
Seq.tOffset = cat(2, zeros(1, Seq.T2Prepare.nEcho+1, nTau1), Seq.tOffset);
if ~isemptyfield(Seq, 'CLTime')
  if all(Seq.CLTime(1)==Seq.CLTime(:))
    % uniform command load time at all tReps
    % use the same for preparation CPMG echo train
    Seq.CLTime = cat(2, ones(1, Seq.T2Prepare.nEcho+1, nTau1)*Seq.CLTime(1), Seq.CLTime);
  else
    % try to mimick CLTime of actual CPMG echo train
    Seq.CLTime = cat(2, ...
      Seq.CLTime(1,2:3,:), ...
      ones(1, Seq.T2Prepare.nEcho-1, nTau1)*Seq.CLTime(1,4,1), ...
      Seq.CLTime);
  end
end


%% TX
PulseExcitation.Bandwidth = Seq.tExcitationBW;
PulseExcitation.FlipAngle = HW.tFlip90Def  * (Seq.excitationAngle/90*pi) / ...
  (HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef);
PulseExcitation.MaxNumberOfSegments = 51;  % FIXME: Why?
PulseExcitation.Frequency = Seq.fLarmor;
PulseExcitation.Phase = 0;
PulseExcitation.iDevice = iDevice;
Seq.T2Prepare.pulseExcitation = Seq.T2Prepare.excitationPulse(HW, 0, PulseExcitation);

PulseRefocus.Bandwidth = Seq.tRefocusBW;
PulseRefocus.FlipAngle = HW.tFlip180Def * (Seq.refocusingAngle/180*pi) / ...
  (HW.TX(iDevice).Amp2FlipPiIn1Sec/HW.TX(iDevice).AmpDef);
PulseRefocus.MaxNumberOfSegments = 51;  % FIXME: Why?
PulseRefocus.Frequency = Seq.fLarmor;
PulseRefocus.Phase = 90;
PulseRefocus.iDevice = iDevice;
Seq.T2Prepare.pulseRefocus = Seq.T2Prepare.refocusingPulse(HW, 0, PulseRefocus);

nActTX = max([...
  size(Seq.T2Prepare.pulseExcitation.Start, 1), ...
  size(Seq.T2Prepare.pulseRefocus.Start, 1), ...
  size(TX.Start, 1)]);

% append NaNs to match sizes
if size(Seq.T2Prepare.pulseExcitation.Start, 1) < nActTX
  Seq.T2Prepare.pulseExcitation.Start(end+1:nActTX,:) = NaN;
  Seq.T2Prepare.pulseExcitation.Duration(end+1:nActTX,:) = NaN;
  Seq.T2Prepare.pulseExcitation.Amplitude(end+1:nActTX,:) = NaN;
  Seq.T2Prepare.pulseExcitation.Frequency(end+1:nActTX,:) = NaN;
  Seq.T2Prepare.pulseExcitation.Phase(end+1:nActTX,:) = NaN;
end
if size(Seq.T2Prepare.pulseRefocus.Start, 1) < nActTX
  Seq.T2Prepare.pulseRefocus.Start(end+1:nActTX,:) = NaN;
  Seq.T2Prepare.pulseRefocus.Duration(end+1:nActTX,:) = NaN;
  Seq.T2Prepare.pulseRefocus.Amplitude(end+1:nActTX,:) = NaN;
  Seq.T2Prepare.pulseRefocus.Frequency(end+1:nActTX,:) = NaN;
  Seq.T2Prepare.pulseRefocus.Phase(end+1:nActTX,:) = NaN;
end
if size(TX.Start, 1) < nActTX
  TX.Start(end+1:nActTX,:) = NaN;
  TX.Duration(end+1:nActTX,:) = NaN;
  TX.Amplitude(end+1:nActTX,:) = NaN;
  TX.Frequency(end+1:nActTX,:) = NaN;
  TX.Phase(end+1:nActTX,:) = NaN;
end

% prepend preparation CPMG echo trains to actual pulse sequence
TX.Start = cat(2, ...
  repmat(Seq.T2Prepare.pulseExcitation.Start, [1, 1, nTau1]), ...
  repmat(Seq.T2Prepare.pulseRefocus.Start, [1, Seq.T2Prepare.nEcho, nTau1]), ...
  TX.Start);
TX.Duration = cat(2, ...
  repmat(Seq.T2Prepare.pulseExcitation.Duration, [1, 1, nTau1]), ...
  repmat(Seq.T2Prepare.pulseRefocus.Duration, [1, Seq.T2Prepare.nEcho, nTau1]), ...
  TX.Duration);
TX.Amplitude = cat(2, ...
  repmat(Seq.T2Prepare.pulseExcitation.Amplitude, [1, 1, nTau1]), ...
  repmat(Seq.T2Prepare.pulseRefocus.Amplitude, [1, Seq.T2Prepare.nEcho, nTau1]), ...
  TX.Amplitude);
TX.Frequency = cat(2, ...
  repmat(Seq.T2Prepare.pulseExcitation.Frequency, [1, 1, nTau1]), ...
  repmat(Seq.T2Prepare.pulseRefocus.Frequency, [1, Seq.T2Prepare.nEcho, nTau1]), ...
  TX.Frequency);
TX.Phase = cat(2, ...
  repmat(Seq.T2Prepare.pulseExcitation.Phase+Seq.pulsePrepare.Phase-Seq.T2Prepare.Phase, [1, 1, nTau1]), ...
  repmat(Seq.T2Prepare.pulseRefocus.Phase+Seq.pulsePrepare.Phase-Seq.T2Prepare.Phase, [1, Seq.T2Prepare.nEcho, nTau1]), ...
  TX.Phase);


%% Digital IO
% FIXME: Actually damp coil. For now, just prepend NaNs.
if ~isemptyfield(Seq, 'DigitalIO')
  Seq.DigitalIO.SetTime = cat(2, ...
    NaN(size(Seq.DigitalIO.SetTime,1), Seq.T2Prepare.nEcho+1, nTau1), ...
    Seq.DigitalIO.SetTime);
  Seq.DigitalIO.SetValue = cat(2, ...
    NaN(size(Seq.DigitalIO.SetValue,1), Seq.T2Prepare.nEcho+1, nTau1), ...
    Seq.DigitalIO.SetValue);
end


%% AQ
% prepend NaNs
AQ.Start = cat(2, NaN(size(AQ.Start,1), Seq.T2Prepare.nEcho+1, nTau1), AQ.Start);
AQ.Frequency = cat(2, NaN(size(AQ.Frequency,1), Seq.T2Prepare.nEcho+1, nTau1), AQ.Frequency);
AQ.Phase = cat(2, NaN(size(AQ.Phase,1), Seq.T2Prepare.nEcho+1, nTau1), AQ.Phase);
AQ.fSample = cat(2, NaN(size(AQ.fSample,1), Seq.T2Prepare.nEcho+1, nTau1), AQ.fSample);
AQ.nSamples = cat(2, NaN(size(AQ.nSamples,1), Seq.T2Prepare.nEcho+1, nTau1), AQ.nSamples);
AQ.ResetPhases = cat(2, zeros(1, Seq.T2Prepare.nEcho+1, nTau1), AQ.ResetPhases);


%% Grad

% move rising flank of slice gradient to before preparation CPMG echo train
iGradSlice = HW.Grad(iDevice).Channel2xyzB(HW.Grad(iDevice).Slice.channel);
if numel(Grad(iGradSlice).Time) > 1 && numel(Grad(iGradSlice).Amp) > 1
  risingTime = Grad(iGradSlice).Time(:,1,:);
  risingAmp = Grad(iGradSlice).Amp(:,1,:);
  isRising = cat(1, ...
    true(1, 1, size(Grad(iGradSlice).Time, 3)), ...
    (risingAmp(2:end,1,:) - risingAmp(1:end-1,1,:))>eps());
  risingTime(~isRising) = NaN;
  risingAmp(~isRising) = NaN;
  fallingTime = Grad(iGradSlice).Time(:,1,:);
  fallingAmp = Grad(iGradSlice).Amp(:,1,:);
  fallingTime(isRising) = NaN;
  fallingAmp(isRising) = NaN;
  notFalling = all(isnan(fallingTime), 1);
  fallingTime(:,:,notFalling) = NaN;
  fallingTime(1,:,notFalling) = 0;
  fallingAmp(:,:,notFalling) = NaN;
  fallingAmp(1,:,notFalling) = Grad(iGradSlice).Amp(1,2,notFalling);
  Grad(iGradSlice).Time = cat(2, ...
    risingTime, ...
    repmat(Grad(iGradSlice).Time(:,3,:), 1, Seq.T2Prepare.nEcho, 1), ...
    fallingTime, ...
    Grad(iGradSlice).Time(:,2:end,:));
  Grad(iGradSlice).Amp = cat(2, ...
    risingAmp, ...
    repmat(Grad(iGradSlice).Amp(:,3,:), 1, Seq.T2Prepare.nEcho, 1), ...
    fallingAmp, ...
    Grad(iGradSlice).Amp(:,2:end,:));
  % Remove leading NaN from start of tReps
  [Grad(iGradSlice).Time, idxSort] = sort(Grad(iGradSlice).Time, 1);
  idxSort = bsxfun(@plus, idxSort, ...
    reshape(((1:size(Grad(iGradSlice).Time, 2)*size(Grad(iGradSlice).Time, 3))-1)*size(Grad(iGradSlice).Time, 1), 1, size(Grad(iGradSlice).Time, 2), size(Grad(iGradSlice).Time, 3)));
  Grad(iGradSlice).Amp = Grad(iGradSlice).Amp(idxSort);
  isAllNaNGrad = all(all(isnan(Grad(iGradSlice).Time),2),3);
  Grad(iGradSlice).Time(isAllNaNGrad,:,:) = [];
  Grad(iGradSlice).Amp(isAllNaNGrad,:,:) = [];
end
if numel(Grad(iGradSlice).Repeat) > 1
  Grad(iGradSlice).Repeat = cat(2, ...
    Grad(iGradSlice).Repeat(1,1,:), ...
    repmat(Grad(iGradSlice).Repeat(1,2,:), 1, Seq.T2Prepare.nEcho, 1), ...
    Grad(iGradSlice).Repeat);
end


% FIXME: How to tread other gradient channels? For now, just prepend NaNs.
for iGrad = setdiff(1:numel(Grad), iGradSlice)
  if numel(Grad(iGrad).Time) > 1
    Grad(iGrad).Time = cat(2, ...
      NaN(size(Grad(iGrad).Time,1), Seq.T2Prepare.nEcho+1, size(Grad(iGrad).Time,3)), ...
      Grad(iGrad).Time);
  end
  if numel(Grad(iGrad).Amp) > 1
    Grad(iGrad).Amp = cat(2, ...
      NaN(size(Grad(iGrad).Amp,1), Seq.T2Prepare.nEcho+1, size(Grad(iGrad).Amp,3)), ...
      Grad(iGrad).Amp);
  end
  if numel(Grad(iGrad).Repeat) > 1
    Grad(iGrad).Repeat = cat(2, ...
      [zeros(1, 1, size(Grad(iGrad).Repeat,3)), ones(1, Seq.T2Prepare.nEcho, size(Grad(iGrad).Repeat,3))], ...
      Grad(iGrad).Repeat);
  end
end

end
