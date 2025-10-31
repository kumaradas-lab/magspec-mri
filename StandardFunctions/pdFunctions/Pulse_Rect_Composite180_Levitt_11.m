function [pulseData] = Pulse_Rect_Composite180_Levitt_11(HW, Center, Pulse, varargin)
%% Create a composite 180 degree RF pulse as described in Levitt et al. [11]
%
%   pulseData = Pulse_Rect_Composite180_Levitt_11(HW, Center, Pulse)
% or:
%   pulseData = Pulse_Rect_Composite180_Levitt_11(HW, Center, Bandwidth, FlipAngle, MaxNumberOfSegments,  MaxLength, Frequency, Phase)
% additionally:
%   excitationAngleFactor = Pulse_Rect_Composite180_Levitt_11(HW, 'Amp')
%   bandwidthFactor = Pulse_Rect_Composite180_Levitt_11(HW, 'Time')
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
%             Maximum number of segments of the pulse (default: 51).
%     MaxLength
%             Maximum length of the pulse in seconds (default: Inf).
%     Frequency
%             Frequency of the rf pulse in Hz (default: HW.fLarmor).
%     Phase   "Local" RF phase of the pulse with respect to the overall sequence
%             in degrees (default: 0).
%     Bandwidth
%             Bandwidth of the pulse in Hz (default:
%             max(FactorTime/Pulse.MaxLength, 2e3) where "FactorTime" is the
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
% (C) Copyright 2017-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% composite settings
FlipAngleComposite = [270; 360;  90; 270; 360;  90];
FlipPhaseComposite = [180;   0;  90; 270;  90;   0];

%% check input
if nargin == 2
  % short path (additional syntax)
  if strcmp(Center, 'Amp')
    pulseData = sum(FlipAngleComposite)/180;
  elseif strcmp(Center, 'Time')
    pulseData = sum(FlipAngleComposite)/180;
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

%% composite pulse
Pulse.FlipAngleComposite = FlipAngleComposite;
Pulse.FlipPhaseComposite = FlipPhaseComposite;
pulseData = Pulse_Rect_Composite180(HW, Center, Pulse);

end
