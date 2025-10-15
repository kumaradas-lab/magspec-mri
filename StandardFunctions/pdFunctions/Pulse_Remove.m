function pulseData = Pulse_Remove(~, Center, ~, varargin)
%% Remove a pulse from pulse program
%
%   pulseData = Pulse_Remove(HW, Center, Pulse)
% or:
%   pulseData = Pulse_Remove(HW, Center, Bandwidth, FlipAngle, MaxNumberOfSegments,  MaxLength, Frequency, Phase)
% additionally:
%   excitationAngleFactor = Pulse_Remove(HW, 'Amp')
%   bandwidthFactor = Pulse_Remove(HW, 'Time')
%
% Independent on any input, this pulse function removes the respective pulse
% from the pulse program completely.
%
% ATTENTION:  Not all sequences support removing a pulse completely using this
%             pulse shape function!
%
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
% (C) Copyright 2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

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

%% dummy pulse
pulseData.Start = [];
pulseData.Amplitude = [];
pulseData.Duration = [];
pulseData.Frequency = [];
pulseData.Phase = [];

end
