function pulseData = Pulse_Rect_Slice(HW, Center, Pulse, varargin)
%% Pulse shape function for a rectangular RF pulse with a matched bandwidth
%
%   pulseData = Pulse_Rect_Slice(HW, Center, Pulse)
% or:
%   pulseData = Pulse_Rect_Slice(HW, Center, Bandwidth, FlipAngle, MaxNumberOfSegments, MaxLength, Frequency, Phase)
% additionally:
%   excitationAngleFactor = Pulse_Rect_Slice(HW, 'Amp')
%   bandwidthFactor = Pulse_Rect_Slice(HW, 'Time')
%
% This function can be used to create rectangular pulses. The bandwidth of these
% pulses is adapted such that the slice thickness corresponds to the FWHM of the
% intensity sinc^2 profile.
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
% (C) Copyright 2018-2025 Pure Devices GmbH, Wuerzburg, Germany
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
bwFactor = 0.442946; % solve sinc(pi*x) = 1/sqrt(2) for x -> x=0.442946
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
% Use gamma that better matches the frequency of the pulse
% FIXME: Could this be an issue with (very) off-center slice pulses?
if abs(Pulse.Frequency - HW.fLarmorX) < abs(Pulse.Frequency - HW.fLarmor)
  tFlipPi = pi/HW.GammaX / HW.TX(Pulse.iDevice).AmpDef;
else
  tFlipPi = HW.TX(Pulse.iDevice).Amp2FlipPiIn1Sec / HW.TX(Pulse.iDevice).AmpDef;
end

BlockLength = 1/Pulse.Bandwidth*bwFactor;

gain = HW.TX(Pulse.iDevice).AmpDef * 2*tFlipPi * (Pulse.FlipAngle/Pulse.FlipAngleFullTurn) / BlockLength;

if Pulse.MaxLength + 1/HW.TX(Pulse.iDevice).fSample < BlockLength
  error('PD:Pulse_Rect_Slice:MaxLengthTooShort', ...
    'MaxLength of rf pulse is %.3f %cs too short.', ...
    (BlockLength - Pulse.MaxLength)*1e6, char(181));
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
