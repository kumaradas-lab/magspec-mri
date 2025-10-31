function pulseData = Pulse_Rect_SpinLock_T1rho(HW, Center, Pulse, varargin)
%% Create two rect RF pulse connected by a spin-locking pulse
%
%   pulseData = Pulse_Rect_SpinLock_T1rho(HW, Center, Pulse)
% or:
%   pulseData = Pulse_Rect_SpinLock_T1rho(HW, Center, Bandwidth, FlipAngle, MaxNumberOfSegments,  maxLength, Frequency, Phase)
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
%             *center* of the second excitation pulse. (Default: 1e-3)
%     AmplitudeSpinLock
%             Amplitude of the spin locking pulse in Tesla.
%             (Default: HW.TX.AmpDef/10)
%     PhaseSpinLock
%             Phase offset of the spin-locking pulse with respect to the "main"
%             excitation pulse in degrees. (Default: 90)
%     cycleB1
%             Invert the B1 pulse for the first half of the spin-locking pulse
%             to reduce the influence of B1 inhomogeneities. [1] (Default: true)
%     withInvert
%             Invert spin system at the center of the spin-locking pulse to
%             reduce influence of B0 inhomogeneities. [1] (Default: true)
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
%
% ------------------------------------------------------------------------------
% (C) Copyright 2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------


%% check input
if nargin == 2
  % short path (additional syntax)
  if strcmp(Center, 'Amp')
    pulseData = 1;
  elseif strcmp(Center, 'Time')
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
  Pulse.AmplitudeSpinLock = HW.TX(Pulse.iDevice).AmpDef/10;
end
if isemptyfield(Pulse, 'PhaseSpinLock'), Pulse.PhaseSpinLock = -90; end
if isemptyfield(Pulse, 'cycleB1'), Pulse.cycleB1 = 1; end
Pulse.cycleB1 = double(Pulse.cycleB1~=0);
if isemptyfield(Pulse, 'withInvert'), Pulse.withInvert = true; end

%% rect pulses for a solid echo
tFlipPi = HW.TX(Pulse.iDevice).Amp2FlipPiIn1Sec / HW.TX(Pulse.iDevice).AmpDef;

BlockLength = 1/Pulse.Bandwidth*0.999;

gain = HW.TX(Pulse.iDevice).AmpDef * 2*tFlipPi * ...
  (Pulse.FlipAngle/Pulse.FlipAngleFullTurn) / (BlockLength/0.998);

if Pulse.MaxLength < BlockLength
  error('MaxLength of rf pulse is too short.');
end

if Pulse.withInvert
  pulseData.Start = ...
    [repmat([-BlockLength/2; BlockLength/2], 1, size(Pulse.DurationSpinLock, 2)); ...
    Pulse.DurationSpinLock/2-BlockLength/2; Pulse.DurationSpinLock/2+BlockLength/2; ...
    Pulse.DurationSpinLock-BlockLength/2] + Center;
  pulseData.Duration = ...
    [repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2)); ...
    (Pulse.DurationSpinLock-2*BlockLength)/2; ...
    repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2)); ...
    (Pulse.DurationSpinLock-2*BlockLength)/2; ...
    repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2))];
  pulseData.Amplitude = [gain; Pulse.AmplitudeSpinLock; 2*gain; Pulse.AmplitudeSpinLock; 0*gain];
  pulseData.Frequency = zeros(size(pulseData.Start)) + Pulse.Frequency;
  % FIXME: Phase-cycle inversion pulse?
  pulseData.Phase = [0; (2*Pulse.cycleB1-1)*Pulse.PhaseSpinLock; Pulse.PhaseSpinLock; ...
    -Pulse.PhaseSpinLock; 180] + Pulse.Phase;
else
  pulseData.Start = ...
    [repmat([-BlockLength/2; BlockLength/2], 1, size(Pulse.DurationSpinLock, 2)); ...
    Pulse.DurationSpinLock/2; Pulse.DurationSpinLock-BlockLength/2] + Center;
  pulseData.Duration = ...
    [repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2)); ...
    (Pulse.DurationSpinLock-BlockLength)/2; ...
    (Pulse.DurationSpinLock-BlockLength)/2; ...
    repmat(BlockLength, 1, size(Pulse.DurationSpinLock, 2))];
  pulseData.Amplitude = [gain; Pulse.AmplitudeSpinLock; Pulse.AmplitudeSpinLock; 0*gain];
  pulseData.Frequency = zeros(size(pulseData.Start)) + Pulse.Frequency;
  pulseData.Phase = [0; (2*Pulse.cycleB1-1)*Pulse.PhaseSpinLock; -Pulse.PhaseSpinLock; 180] + Pulse.Phase;
end

end
