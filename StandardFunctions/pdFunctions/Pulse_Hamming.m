function [pulseData] = Pulse_Hamming(HW, Center, Pulse, varargin)
%% Create a "rect" pulse with Hamming window
%
%   pulseData = Pulse_Hamming(HW, Center, Pulse)
% or:
%   pulseData = Pulse_Hamming(HW, Center, Bandwidth, FlipAngle, MaxNumberOfSegments,  maxLength, Frequency, Phase)
% additionally:
%   excitationAngleFactor = Pulse_Hamming(HW, 'Amp')
%   bandwidthFactor = Pulse_Hamming(HW, 'Time')
%
% This function is used for creating shaped rf pulses which correspond to a rect
% pulse filtered by a Hamming window (see Hamming.m).
% Simplifying, the Hamming window reduces the amplitude of the side bands that
% are excited by the pulse.
%
% INPUT:
%   HW      HW structure
%   Center  The center of the pulse in seconds (tRep).
%   Pulse   A structure with the following fields (if omitted or empty, default
%           values are used):
%     FlipAngle
%             The total (effective) flip angle of the pulse in radians (or the
%             units defined by FlipAngleFullTurn). It is used to set an
%             appropriate pulse amplitude (default: pi).
%     FlipAngleFullTurn
%             Value that defines a full turn in FlipAngle units (e.g. 360 for
%             degrees, or 2*pi for radians, default: 2*pi).
%     MaxNumberOfSegments
%             Maximum number of segments of the pulse (default: 51).
%     MaxLength
%             Maximum length of the pulse in seconds (default: Inf).
%     Frequency
%             Frequency of the rf pulse in Hz (default: HW.fLarmor).
%     Phase   "Local" RF phase of the pulse with respect to the overall sequence
%             in degrees (default: 0).
%     Bandwidth
%             Bandwidth of the pulse in Hz
%             (default: max(1/Pulse.MaxLength, 2e3) )
%
% OUTPUT:
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
% The additional syntax is used to return the amplitude and bandwidth factors.
% The "excitationAngleFactor" is the factor that must be applied to the
% amplitude of this pulse to have the same excitation angle as a rect pulse of
% the same length. The "bandwidthFactor" is the factor that must be applied to
% the duration of the pulse to have the same bandwidth (FWHM) as a rect pulse.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

%% check input
if nargin == 2
  % short path (additional syntax)
  if strcmp(Center, 'Amp')
    pulseData = 2;  %1/0.54348;
  elseif strcmp(Center, 'Time')
    pulseData = 2;  %1/0.54348;
  else
    pulseData = NaN;
  end
  return;
end

%% Convert from syntax (2) to syntax (1)
if ~isstruct(Pulse), Pulse = struct('Bandwidth', Pulse); end
if nargin > 3, Pulse.FlipAngle = varargin{1}; end
if nargin > 4, Pulse.MaxNumberOfSegments = varargin{2}; end
if nargin > 5, Pulse.MaxLength = varargin{3}; end
if nargin > 6, Pulse.Frequency = varargin{4}; end
if nargin > 7, Pulse.Phase = varargin{5}; end

%% default values
Pulse = set_EmptyField(Pulse, 'FlipAngle', pi);
Pulse = set_EmptyField(Pulse, 'FlipAngleFullTurn', 2*pi);
Pulse = set_EmptyField(Pulse, 'MaxNumberOfSegments', 51);
Pulse = set_EmptyField(Pulse, 'MaxLength', Inf);
Pulse = set_EmptyField(Pulse, 'Frequency', HW.fLarmor);
Pulse = set_EmptyField(Pulse, 'Phase', 0);
Pulse = set_EmptyField(Pulse, 'Bandwidth', max(1/Pulse.MaxLength, 2e3));  % FIXME: Is this a reasonable default?
Pulse = set_EmptyField(Pulse, 'iDevice', 1);

%% Hamming filtered pulse
numberOfZeroCrossings = 2;  % Number ot times the shaped pulse crosses zero [2,4,6,8,...]

pulseDuration = numberOfZeroCrossings / Pulse.Bandwidth;  % duration of the whole pulse
if Pulse.MaxLength + 1/HW.TX(Pulse.iDevice).fSample < pulseDuration
  error('PD:Pulse_Hamming:MaxLengthTooShort', 'MaxLength is too short');
end

% reduce the number of segments if the pulse is too short
if Pulse.MaxNumberOfSegments > pulseDuration/(2/HW.TX(Pulse.iDevice).fSample)
  % one segment should be at least 2 HF DAC samples
  numberOfSegments = floor(pulseDuration/(2/HW.TX(Pulse.iDevice).fSample));
else
  % if the pulse is long enough, use maxNumbersOfSegments
  numberOfSegments = Pulse.MaxNumberOfSegments;
end

% timeline rounded to DAC sample times
tStart = round((linspace(-pulseDuration/2, ...
                         pulseDuration/2, ...
                         numberOfSegments+1).' + Center) ...
               * HW.TX(Pulse.iDevice).fSample) / HW.TX(Pulse.iDevice).fSample;

% calculate the start time such that the tShape time is centered on the rf pulse
pulseData.Start = tStart(1:end-1);
% duration of pulse segments
pulseData.Duration = diff(tStart);
% frequency of pulse segments
pulseData.Frequency = zeros(numberOfSegments, 1) + Pulse.Frequency;
HamWin = Hamming(numberOfSegments*2+1);
% normalized amplitude at tShape
B1Shape = HamWin(2:2:end);
% amplitude (in Tesla) of the B1+ field in the coil
% Use gamma that better matches the frequency of the pulse
% FIXME: Could this be an issue with (very) off-center slice pulses?
if abs(Pulse.Frequency - HW.fLarmorX) < abs(Pulse.Frequency - HW.fLarmor)
  pulseData.Amplitude = abs(B1Shape) * ...
    (2*pi/HW.GammaX * Pulse.FlipAngle/Pulse.FlipAngleFullTurn) / ...
    sum(pulseData.Duration.*B1Shape, 1);
else
  pulseData.Amplitude = abs(B1Shape) * ...
    (2*HW.TX(Pulse.iDevice).Amp2FlipPiIn1Sec * Pulse.FlipAngle/Pulse.FlipAngleFullTurn) / ...
    sum(pulseData.Duration.*B1Shape, 1);
end
% Phase of pulse segments (no negative Amplitude is allowed, so you have to add
% 180 degrees to the phase to get an equivalent)
pulseData.Phase = angle(B1Shape)/pi*180 + 0 + Pulse.Phase;

end
