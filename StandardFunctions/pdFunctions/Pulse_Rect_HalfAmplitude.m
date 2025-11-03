function [pulseData] = Pulse_Rect_HalfAmplitude(HW, Center, Pulse, varargin)
%% Create a rectangular RF pulse with TX.Amplitude = HW.TX.AmpDef / 2
%
%   pulseData = Pulse_Rect_HalfAmplitude(HW, Center, Pulse)
% or:
%   pulseData = Pulse_Rect_HalfAmplitude(HW, Center, Bandwidth, FlipAngle, MaxNumberOfSegments,  maxLength, Frequency, Phase)
% additionally:
%   excitationAngleFactor = Pulse_Rect_HalfAmplitude(HW, 'Amp')
%   bandwidthFactor = Pulse_Rect_HalfAmplitude(HW, 'Time')
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
% (C) Copyright 2012-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

%% bandwidth factor
% x = 2 * 6/5/0.86;
% x = 2 * 12/9/0.86;
% x = 2 * 12/9/0.86;
% x = 2*6/4/0.86;
x = 2;

%% check input
if nargin == 2
  % short path (additional syntax)
  if strcmp(Center, 'Amp')
    pulseData =x ;
  elseif strcmp(Center, 'Time')
    pulseData = 1/x;
  else
    pulseData = NaN;
  end
  return;
end

tFlipPi = HW.TX(Pulse.iDevice).Amp2FlipPiIn1Sec / HW.TX(Pulse.iDevice).AmpDef;

BlockLength=1/x/Pulse.Bandwidth*0.999;

gain = HW.TX(Pulse.iDevice).AmpDef * 2*tFlipPi * ...
  (Pulse.FlipAngle/Pulse.FlipAngleFullTurn) / (BlockLength/0.998);

if Pulse.MaxLength < BlockLength
  error('MaxLength of rf pulse is too short.');
end

if Pulse.MaxNumberOfSegments < 1
  error('MaxNumberOfSegments must be at least 1.');
end

pulseData.Start = -BlockLength/2 + Center;
pulseData.Amplitude = gain;
pulseData.Duration = BlockLength;
pulseData.Frequency = Pulse.Frequency;
pulseData.Phase = Pulse.Phase;

end
