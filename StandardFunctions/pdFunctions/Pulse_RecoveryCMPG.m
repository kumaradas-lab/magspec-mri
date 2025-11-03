function pulseData = Pulse_RecoveryCMPG(HW, Center, Pulse, varargin)
%% Pulse shape function for rectangular pulse with the sequence_RecoveryCPMG
%
%   pulseData = Pulse_RecoveryCMPG(HW, Center, Pulse)
% or:
%   pulseData = Pulse_RecoveryCMPG(HW, Center, Bandwidth, FlipAngle, MaxNumberOfSegments,  MaxLength, Frequency, Phase)
% additionally:
%   excitationAngleFactor = Pulse_Rect_Slice(HW, 'Amp')
%   bandwidthFactor = Pulse_RecoveryCMPG(HW, 'Time')
%
% This function can be used to create rectangular pulses. It is solely to be
% used with sequence_RecoveryCPMG where the following applies (as verified by
% measurements):
% For a profile of bandwidth BW (measured between the 1/sqrt(2) of the max.
% amplitude values in Hz), use the following combination of pulse durations and
% acquisition window length:
% p180 = 1/BW / 2;
% p90 = p180 * 7/4;
% tAQ = p180 * 3/4;
%
% INPUT:
%   HW        HW structure
%   Center    The center of the pulse in seconds (tRep).
%   Pulse     A structure with the following fields (if omitted or empty,
%             default values are used):
%     FlipAngle
%             The total (effective) flip angle of the pulse in radians (or the
%             units defined by FlipAngleFullTurn). It is used to set an
%             appropriate pulse amplitude (default: pi/2).
%     FlipAngleFullTurn
%             Value that defines a full turn in FlipAngle units (e.g. 360 for
%             degrees, or 2*pi for radians, default: 2*pi).
%     MaxNumberOfSegments
%             Maximum number of segments of the pulse. Ignored for rectangular
%             pulse which always consists of one segment (default: Inf).
%     MaxLength
%             Maximum length of the pulse in seconds (default: Inf).
%     Frequency
%             Frequency of the rf pulse in Hz (default: HW.fLarmor).
%     Phase   "Local" RF phase of the pulse with respect to the overall sequence
%             in degrees (default: 0).
%     Bandwidth
%             Bandwidth of the pulse in Hz (default:
%             max(1/Pulse.MaxLength, 2e3) where "FactorTime" is the
%             multiplier compared to a simple rectangular pulse with the same
%             flip angle.
%
% OUTPUT:
%   pulseData A structure with the following fields:
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
% The additional syntax is used to return the amplitude and bandwidth factors.
% The "excitationAngleFactor" is the factor that must be applied to the
% amplitude of this pulse to have the same excitation angle as a rect pulse of
% the same length. The "bandwidthFactor" is the factor that must be applied to
% the duration of the pulse to have the same bandwidth (FWHM) as a rect pulse.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2019-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% Convert from syntax (2) to syntax (1)
if nargin < 3, Pulse = []; end
if ~isstruct(Pulse), Pulse = struct('Bandwidth', Pulse); end
if nargin > 3, Pulse.FlipAngle = varargin{1}; end
if nargin > 4, Pulse.MaxNumberOfSegments = varargin{2}; end
if nargin > 5, Pulse.MaxLength = varargin{3}; end
if nargin > 6, Pulse.Frequency = varargin{4}; end
if nargin > 7, Pulse.Phase = varargin{5}; end

%% default values
bwFactor = 0.5*1.28;  % This factor is derived from experiments with the magspec.
% FIXME: Maybe we'd need a different factor for the pulse program with 2*tAQ.
Pulse = set_EmptyField(Pulse, 'FlipAngle', pi/2);
Pulse = set_EmptyField(Pulse, 'FlipAngleFullTurn', 2*pi);
Pulse = set_EmptyField(Pulse, 'MaxNumberOfSegments', Inf);
Pulse = set_EmptyField(Pulse, 'MaxLength', Inf);
Pulse = set_EmptyField(Pulse, 'Frequency', HW.fLarmor);
Pulse = set_EmptyField(Pulse, 'Phase', 0);
Pulse = set_EmptyField(Pulse, 'Bandwidth', max(1/Pulse.MaxLength/bwFactor, 2e3)); % FIXME: Is this a sensible default?
Pulse = set_EmptyField(Pulse, 'iDevice', 1);

%% short path
if nargin == 2
  if ischar(Center)
    if strcmp(Center, 'Amp')
      pulseData = 1;
    elseif strcmp(Center, 'Time')
      pulseData = bwFactor;
    else
      pulseData = NaN;
    end
  else
    pulseData = NaN;
  end
  return;
end

%% "real" path
tFlipPi = HW.TX(Pulse.iDevice).Amp2FlipPiIn1Sec/HW.TX(Pulse.iDevice).AmpDef;

BlockLength = 1/Pulse.Bandwidth*bwFactor;

gain = HW.TX(Pulse.iDevice).AmpDef * 2*tFlipPi * (Pulse.FlipAngle/Pulse.FlipAngleFullTurn) / BlockLength;

if Pulse.MaxLength < BlockLength;
  error('maxLength of HF Pulse to short')
end

if Pulse.MaxNumberOfSegments < 1
  error('Pulse.MaxNumberOfSegments must be >= 1');
end

pulseData.Start = -BlockLength/2+Center;
pulseData.Amplitude = gain;
pulseData.Duration = BlockLength;
pulseData.Frequency = Pulse.Frequency;
pulseData.Phase = Pulse.Phase;

end
