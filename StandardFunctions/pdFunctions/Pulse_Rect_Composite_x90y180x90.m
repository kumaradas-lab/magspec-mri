function [pulseData] = Pulse_Rect_Composite_x90y180x90(HW, Center, Pulse, varargin)
%% Create a composite 180 degree inversion RF pulse
%
%   pulseData = Pulse_Rect_Composite_x90y180x90(HW, Center, Pulse)
% or:
%   pulseData = Pulse_Rect_Composite_x90y180x90(HW, Center, Bandwidth, FlipAngle, MaxNumberOfSegments,  maxLength, Frequency, Phase)
% additionally:
%   excitationAngleFactor = Pulse_Rect_Composite_x90y180x90(HW, 'Amp')
%   bandwidthFactor = Pulse_Rect_Composite_x90y180x90(HW, 'Time')
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
% (C) Copyright 2016-2020 Pure Devices GmbH, Wuerzburg, Germany
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

%% composite pulse
tFlipPi = HW.TX(Pulse.iDevice).Amp2FlipPiIn1Sec / HW.TX(Pulse.iDevice).AmpDef;

BlockLength = 1/Pulse.Bandwidth * 0.999;

gain = HW.TX(Pulse.iDevice).AmpDef * 2*tFlipPi * ...
  (Pulse.FlipAngle/Pulse.FlipAngleFullTurn) / (BlockLength/0.998);

if Pulse.MaxLength < BlockLength
  error('PD:Pulse_Rect_Composite_x90y180x90:MaxLengthTooShort', 'MaxLength is too short');
end

if Pulse.MaxNumberOfSegments < 3
  error('MaxNumberOfSegments must be at least 3.');
end

pulseData.Start = [-BlockLength; -BlockLength/2; BlockLength/2] + Center;
pulseData.Amplitude = [gain; gain; gain];
pulseData.Duration = [BlockLength/2; BlockLength; BlockLength/2];
pulseData.Frequency= [Pulse.Frequency; Pulse.Frequency; Pulse.Frequency];
pulseData.Phase= 45 + [0;90;0] + Pulse.Phase - 90;

end
