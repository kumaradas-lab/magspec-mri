function pulseData = Pulse_Rect_SpinLock_T1rho(HW, Center, Pulse, varargin)
%% Create two rect RF pulse connected by a spin-locking pulse
%
%   pulseData = Pulse_Rect_SpinLock_T1rho(HW, Center, Pulse)
% or:
%   pulseData = Pulse_Rect_SpinLock_T1rho(HW, Center, Bandwidth, FlipAngle, MaxNumberOfSegments, MaxLength, Frequency, Phase)
% additionally:
%   excitationAngleFactor = Pulse_Rect_SpinLock_T1rho(HW, 'Amp')
%   bandwidthFactor = Pulse_Rect_SpinLock_T1rho(HW, 'Time')
%
% A spin-locking pulse is a continuous "rf pulse" that can be used to determine
% relaxation parameters in the rotating frame of reference of the B1 pulse.
% The spin-locking pulse is immediately adjacent to the end of the first
% excitation pulse and the start of the second excitation pulse.
%
%
% INPUT:
%
%   HW      HW structure
%   Center  The center of the first excitation pulse in seconds (tRep).
%   Pulse   A structure with the following fields (if omitted or empty, default
%           values are used):
%     FlipAngle
%             The total (effective) flip angle of the pulse in radians (or the
%             units defined by FlipAngleFullTurn). It is used to set an
%             appropriate pulse amplitude (default: pi/2).
%     FlipAngleFullTurn
%             Value that defines a full turn in FlipAngle units (e.g. 360 for
%             degrees, or 2*pi for radians, default: 2*pi).
%     MaxNumberOfSegments (unused)
%             Maximum number of segments of the pulse (default: 1).
%     MaxLength
%             Maximum length of the pulse in seconds (default: Inf).
%     Frequency
%             Frequency of the rf pulse in Hz (default: HW.fLarmor).
%     Phase
%             "Local" RF phase of the pulse with respect to the overall sequence
%             in degrees (default: 0).
%     Bandwidth
%             Bandwidth of the pulse in Hz
%             (default: max(1/Pulse.MaxLength, 2e3) )
%     DurationSpinLock
%             "Duration" of the spin-locking pulse in seconds. Instead of the
%             actual duration of the spin-locking pulse, this is the time
%             between the *center* of the first excitation pulse and the
%             *center* of the (omitted) second excitation pulse. (Default: 1e-3)
%     AmplitudeSpinLock
%             Amplitude of the spin locking pulse in Hertz.
%             (Default: HW.TX.AmpDef/10*2*pi/HW.GammaDef)
%     PhaseSpinLock
%             Phase offset of the spin-locking pulse with respect to the "main"
%             excitation pulse in degrees. (Default: 90)
%     cycleB1
%             Invert the B1 pulse for the first half of the spin-locking pulse
%             to reduce the influence of B1 inhomogeneities. [1] (Default: true)
%     withInvert
%             Method for spin system inversion:
%             0: rotary echo spin-lock, i.e., invert the phase of the spin-lock
%                after half its duration.
%             1: composite spin-lock, i.e., invert the spin system at the center
%                of the spin-locking pulse to reduce influence of B0
%                inhomogeneities. [1]
%             2: balanced spin lock, i.e., invert the spin system at one quarter
%                and three quarters of its duration to reduce influence of B0
%                and B1 homogeneities. [2]
%             (Default: 1)
%     withFinalFlip
%             Boolean value that selects whether the magnetization is flipped
%             back to the z-direction after the spin-lock pulse.
%             (Default: false)
%
%
% OUTPUT:
%
%   pulseData
%          A structure with the following fields:
%     Start   Column vector with the start times of each component/block in
%             seconds (tRep).
%     Amplitude
%             Column vector with the amplitudes of each component/block in
%             Tesla.
%     Duration
%             Column vector with the durations of each component/block in
%             seconds.
%     Frequency
%             Column vector with the frequencies of each component/block in Hz.
%     Phase
%             Column vector with the phases in degrees of each component/block
%             with respect to the overall sequence.
%
%
% The additional syntax is used to return the amplitude and bandwidth factors.
% The "excitationAngleFactor" is the factor that must be applied to the
% amplitude of this pulse to have the same excitation angle as a rect pulse of
% the same length. The "bandwidthFactor" is the factor that must be applied to
% the duration of the pulse to have the same bandwidth (FWHM) as a rect pulse.
%
%
% [1]: https://qims.amegroups.com/article/view/8702/9369
% [2]: https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28585
%
% ------------------------------------------------------------------------------
% (C) Copyright 2021-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% check input
if nargin == 2
  % short path (additional syntax)
  if strcmp(Center, 'Amp')
    % FIXME: Is there a way to check if the inversion pulse will be used?
    % if isemptyfield(Pulse, 'withInvert') || Pulse.withInvert > 0
      % The inversion pulse(s) have twice the amplitude of the excitation pulse.
      pulseData = 2;
    % else
    %   pulseData = 1;
    % end
  elseif strcmp(Center, 'Time')
    % This only makes sense towards before the rf pulse.
    pulseData = 1;
  else
    pulseData = NaN;
  end
  return;
end


%% Convert from syntax (2) to syntax (1)
if nargin < 3, Pulse = []; end
if ~isstruct(Pulse), Pulse = struct('Bandwidth', Pulse); end
if nargin > 3, Pulse.FlipAngle = varargin{1}; end
if nargin > 4, Pulse.MaxNumberOfSegments = varargin{2}; end
if nargin > 5, Pulse.MaxLength = varargin{3}; end
if nargin > 6, Pulse.Frequency = varargin{4}; end
if nargin > 7, Pulse.Phase = varargin{5}; end


%% default values
Pulse = set_EmptyField(Pulse, 'FlipAngle', pi/2);
Pulse = set_EmptyField(Pulse, 'FlipAngleFullTurn', 2*pi);
Pulse = set_EmptyField(Pulse, 'MaxNumberOfSegments', 1);
Pulse = set_EmptyField(Pulse, 'MaxLength', Inf);
Pulse = set_EmptyField(Pulse, 'Frequency', HW.fLarmor);
Pulse = set_EmptyField(Pulse, 'Phase', 0);
Pulse = set_EmptyField(Pulse, 'Bandwidth', max(1/Pulse.MaxLength, 2e3));  % FIXME: Is this a sensible default?
Pulse = set_EmptyField(Pulse, 'iDevice', 1);
% pulse shape specific settings
if isemptyfield(Pulse, 'DurationSpinLock'), Pulse.DurationSpinLock = 1e-3; end
if isemptyfield(Pulse, 'AmplitudeSpinLock')
  Pulse.AmplitudeSpinLock = HW.TX(Pulse.iDevice).AmpDef/10*HW.GammaDef/(2*pi);
end
if isemptyfield(Pulse, 'PhaseSpinLock'), Pulse.PhaseSpinLock = -90; end
if isemptyfield(Pulse, 'cycleB1'), Pulse.cycleB1 = 1; end
Pulse.cycleB1 = double(Pulse.cycleB1~=0);
if isemptyfield(Pulse, 'withInvert'), Pulse.withInvert = 1; end
if isemptyfield(Pulse, 'withFinalFlip'), Pulse.withFinalFlip = 0; end


%% rect pulse with spin lock
% Use gamma that better matches the frequency of the pulse
% FIXME: Could this be an issue with (very) off-center slice pulses?
if abs(Pulse.Frequency - HW.fLarmorX) < abs(Pulse.Frequency - HW.fLarmor)
  tFlipPi = pi/HW.GammaX / HW.TX(Pulse.iDevice).AmpDef;
else
  tFlipPi = HW.TX(Pulse.iDevice).Amp2FlipPiIn1Sec / HW.TX(Pulse.iDevice).AmpDef;
end

BlockLength = 1 / Pulse.Bandwidth;

gain = HW.TX(Pulse.iDevice).AmpDef * 2*tFlipPi * ...
  (Pulse.FlipAngle/Pulse.FlipAngleFullTurn) ./ BlockLength;

if Pulse.MaxLength + 1/HW.TX(Pulse.iDevice).fSample < BlockLength(1)
  error('PD:Pulse_Rect_SpinLock_T1rho:MaxLengthTooShort', ...
    'MaxLength (%.3 %cs) of rf pulse is %.3f %cs too short.', ...
    Pulse.MaxLength(1)*1e6, 181, (BlockLength(1) - Pulse.MaxLength)*1e6, 181);
end

if Pulse.withFinalFlip
  finalFlipFactor = 1;
else
  finalFlipFactor = NaN;
end

% convert from Hertz to Tesla
spinlockAmp = Pulse.AmplitudeSpinLock / (HW.GammaDef/(2*pi));

if Pulse.withFinalFlip
  rawSpinLockDuration = Pulse.DurationSpinLock;
else
  rawSpinLockDuration = Pulse.DurationSpinLock+BlockLength/2;
end

switch Pulse.withInvert
  case 0
    pulseData.Start = ...
      [repmat([-BlockLength/2; BlockLength/2], 1, size(Pulse.DurationSpinLock, 2)); ...
      rawSpinLockDuration/2; rawSpinLockDuration-BlockLength/2] + Center;
    pulseData.Duration = ...
      [repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2)); ...
      (rawSpinLockDuration-BlockLength)/2; ...
      (rawSpinLockDuration-BlockLength)/2; ...
      repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2))];
    pulseData.Amplitude = [gain; spinlockAmp; spinlockAmp; finalFlipFactor*gain];
    pulseData.Frequency = zeros(size(pulseData.Start)) + Pulse.Frequency;
    pulseData.Phase = [0; (2*Pulse.cycleB1-1)*Pulse.PhaseSpinLock; -Pulse.PhaseSpinLock; 180] + Pulse.Phase;

  case 1
    pulseData.Start = ...
      [repmat([-BlockLength/2; BlockLength/2], 1, size(Pulse.DurationSpinLock, 2)); ...
      rawSpinLockDuration/2-BlockLength/2;
      rawSpinLockDuration/2+BlockLength/2; ...
      rawSpinLockDuration-BlockLength/2] + Center;
    pulseData.Duration = ...
      [repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2)); ...
      (rawSpinLockDuration-2*BlockLength)/2; ...
      repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2)); ...
      (rawSpinLockDuration-2*BlockLength)/2; ...
      repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2))];
    pulseData.Amplitude = [gain; spinlockAmp; 2*gain; spinlockAmp; finalFlipFactor*gain];
    pulseData.Frequency = zeros(size(pulseData.Start)) + Pulse.Frequency;
    % FIXME: Phase-cycle inversion pulse?
    pulseData.Phase = [0; (2*Pulse.cycleB1-1)*Pulse.PhaseSpinLock; Pulse.PhaseSpinLock; ...
      -Pulse.PhaseSpinLock; 180] + Pulse.Phase;

  case 2
    pulseData.Start = ...
      [repmat([-BlockLength/2; BlockLength/2], 1, size(Pulse.DurationSpinLock, 2)); ...
      rawSpinLockDuration/4-BlockLength/2;
      rawSpinLockDuration/4+BlockLength/2; ...
      3*rawSpinLockDuration/4-BlockLength/2;
      3*rawSpinLockDuration/4+BlockLength/2; ...
      rawSpinLockDuration-BlockLength/2] + Center;
    pulseData.Duration = ...
      [repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2)); ...
      rawSpinLockDuration/4-BlockLength; ...
      repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2)); ...
      rawSpinLockDuration/2-BlockLength; ...
      repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2)); ...
      rawSpinLockDuration/4-BlockLength;
      repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2))];
    pulseData.Amplitude = [gain; spinlockAmp; 2*gain; spinlockAmp; 2*gain; spinlockAmp; finalFlipFactor*gain];
    pulseData.Frequency = zeros(size(pulseData.Start)) + Pulse.Frequency;
    % FIXME: Phase-cycle inversion pulses?
    pulseData.Phase = [0; ...  % excitation
      (2*Pulse.cycleB1-1)*Pulse.PhaseSpinLock; ...  % spin-lock
      Pulse.PhaseSpinLock; ...  % inversion
      -Pulse.PhaseSpinLock; ...  % (minus) spin-lock
      -Pulse.PhaseSpinLock; ...  % inversion
      (2*Pulse.cycleB1-1)*Pulse.PhaseSpinLock; ...  % spin-lock
      180] ...  % excitation
      + Pulse.Phase;

end

pulseData.Start(isnan(pulseData.Amplitude),:) = NaN;
pulseData.Frequency(isnan(pulseData.Amplitude),:) = NaN;
pulseData.Phase(isnan(pulseData.Amplitude(:,1)),:) = NaN;

allNaNRows = all(isnan(pulseData.Amplitude), 2);
pulseData.Start(allNaNRows,:) = [];
pulseData.Duration(allNaNRows,:) = [];
pulseData.Frequency(allNaNRows,:) = [];
pulseData.Phase(allNaNRows,:) = [];
pulseData.Amplitude(allNaNRows,:) = [];

end
