function pulseData = Pulse_Rect_Composite180(HW, Center, Pulse, varargin)
%% Create a composite 180 degree (inversion) RF pulse
%
%   pulseData = Pulse_Rect_Composite180(HW, Center, Pulse)
% or:
%   pulseData = Pulse_Rect_Composite180(HW, Center, Bandwidth, FlipAngle, MaxNumberOfSegments, MaxLength, Frequency, Phase)
% additionally:
%   excitationAngleFactor = Pulse_Rect_Composite180(HW, 'Amp')
%   bandwidthFactor = Pulse_Rect_Composite180(HW, 'Time')
%
% This function can be used to create inversion pulses that are composed of
% several arbitrarily definable rectangular pulses of equal amplitude.  The
% first calling syntax above uses a structure "Pulse" to pass parameters to the
% function.  The second syntax (deprecated) can be used to pass parameters as
% separate arguments.
%
% INPUT:
%   HW        HW structure
%   Center    The center of the pulse in seconds (tRep).
%   Pulse     A structure with the following fields (if omitted or empty,
%             default values are used):
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
%     FlipAngleComposite
%             Column vector with flip angles in degrees for each component/block
%             of the composite pulse (default: [90; 180; 90]).
%     FlipPhaseComposite
%             Column vector with RF phases in degrees for each component/block
%             of the composite pulse with respect to the "local" phase (default:
%             [0; 90; 0]-90).
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
% amplitude of this pulse to have the same inversion angle as a rect pulse of
% the same length. The "bandwidthFactor" is the factor that must be applied to
% the duration of the pulse to have the same bandwidth (FWHM) as a rect pulse.
%
% ------------------------------------------------------------------------
% (C) Copyright 2012-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------

%% Convert from syntax (2) to syntax (1)
if nargin < 3, Pulse = []; end
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
Pulse = set_EmptyField(Pulse, 'FlipAngleComposite', [ 90;180; 90]);
Pulse = set_EmptyField(Pulse, 'FlipPhaseComposite', [  0; 90;  0]-90);
FactorTime = sum(Pulse.FlipAngleComposite)/180;
Pulse = set_EmptyField(Pulse, 'Bandwidth', max(FactorTime/Pulse.MaxLength, 2e3)); % FIXME: Is this a sensible default?
Pulse = set_EmptyField(Pulse, 'iDevice',  1);

if nargin == 2
  %% short path (additional syntax)
  if strcmp(Center, 'Amp')
    pulseData = sum(Pulse.FlipAngleComposite)/180;
  elseif strcmp(Center, 'Time')
    pulseData = FactorTime;
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

BlockLength = 1/Pulse.Bandwidth;

gain = HW.TX(Pulse.iDevice).AmpDef * 2*tFlipPi * (Pulse.FlipAngle/Pulse.FlipAngleFullTurn) / BlockLength;

if Pulse.MaxLength + 1/HW.TX(Pulse.iDevice).fSample < BlockLength*sum(Pulse.FlipAngleComposite)/180
  error('PD:Pulse_Rect_Composite180:MaxLengthTooShort', ...
    'Pulse_Rect_Composite180: MaxLength of rf pulse is %.3f %cs too short.', ...
    (BlockLength*sum(Pulse.FlipAngleComposite)/180 - Pulse.MaxLength)*1e6, char(181));
end

if Pulse.MaxNumberOfSegments < numel(Pulse.FlipAngleComposite)
  error(['Pulse_Rect_Composite180: MaxNumberOfSegments >= ' num2str(numel(Pulse.FlipAngleComposite))]);
end

pulseData.Start = (cumsum([0; Pulse.FlipAngleComposite(1:end-1)])-sum(Pulse.FlipAngleComposite)/2) / 180*BlockLength+Center;
pulseData.Amplitude = ones(size(Pulse.FlipAngleComposite)) * gain;
pulseData.Duration = Pulse.FlipAngleComposite / 180*BlockLength;
pulseData.Frequency = ones(size(Pulse.FlipAngleComposite)) * Pulse.Frequency;
pulseData.Phase = Pulse.FlipPhaseComposite + Pulse.Phase;

end
