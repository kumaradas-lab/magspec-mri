function [pulseData] = Pulse_RaisedCos(HW, Center, Pulse, varargin)
%% Create a pulse with the shape of a raised cosine
%
%   pulseData = Pulse_RaisedCos(HW, Center, Pulse)
% or:
%   pulseData = Pulse_RaisedCos(HW, Center, Bandwidth, FlipAngle, MaxNumberOfSegments,  maxLength, Frequency, Phase)
% additionally:
%   excitationAngleFactor = Pulse_RaisedCos(HW, 'Amp')
%   bandwidthFactor = Pulse_RaisedCos(HW, 'Time')
%
% This function is used for creating rf pulses with the shape of a raised
% cosine. That is the pulse is shape like the center lobe of a squares cosine
% function.
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
% (C) Copyright 2014-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

%% check input
if nargin == 2
  % short path (additional syntax)
  if strcmp(Center, 'Amp')
    pulseData = 2;
  elseif strcmp(Center, 'Time')
    pulseData = 2;
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

%% raised cosine pulse
numberOfZeroCrossings = 2;  % number of times the cos^2 function crosses zero [2,4,6,8,...]

pulseDuration = numberOfZeroCrossings / Pulse.Bandwidth;  % duration of the whole pulse
if Pulse.MaxLength + 1/HW.TX(Pulse.iDevice).fSample < pulseDuration
  error('MaxLength is too short');
end

% reduce the number of segments if the pulse is too short
if Pulse.MaxNumberOfSegments > pulseDuration/(1/HW.TX(Pulse.iDevice).fSample*2)
  % one segment should be at least 2 HF DAC samples
  numberOfSegments = floor(pulseDuration/(1/HW.TX(Pulse.iDevice).fSample*2));
else
  % if the pulse is long enough use maxNumbersOfSegments
  numberOfSegments = Pulse.MaxNumberOfSegments;
end

% timeline rounded to DAC sample times
tShape = round(linspace(-pulseDuration/2 + pulseDuration/numberOfSegments, ...
                        pulseDuration/2 - pulseDuration/numberOfSegments, ...
                        numberOfSegments).' ...
               * (HW.TX(Pulse.iDevice).fSample/2)) / (HW.TX(Pulse.iDevice).fSample/2);

% calculate the start time such that the tShape time is centered on the rf pulse
pulseData.Start = tShape - [diff(tShape); tShape(end)-tShape(end-1)]/2 + Center;
% duration of pulse segments
pulseData.Duration = [diff(pulseData.Start); pulseData.Start(end)-pulseData.Start(end-1)];
% frequency of pulse segments
pulseData.Frequency = zeros(numberOfSegments,1) + Pulse.Frequency;
% normalized amplitude at tShape
B1Shape = cos(pi/2 * tShape / pulseDuration * numberOfZeroCrossings).^2;
% amplitude (in Tesla) of the B1+ field in the coil
pulseData.Amplitude = abs(B1Shape) * ...
  (2*HW.TX(Pulse.iDevice).Amp2FlipPiIn1Sec * Pulse.FlipAngle/Pulse.FlipAngleFullTurn) / ...
  sum(pulseData.Duration.*B1Shape, 1);
% Phase of pulse segments (no negative Amplitude is allowed, so you have to add
% 180 degrees to the phase to get an equivalent)
pulseData.Phase = angle(B1Shape)/pi*180 + 0 + Pulse.Phase;

end
