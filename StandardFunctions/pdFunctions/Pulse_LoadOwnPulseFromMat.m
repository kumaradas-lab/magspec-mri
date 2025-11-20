function [TXShape] = Pulse_LoadOwnPulseFromMat(HW, Center, BW, FlipAngle, maxpulseCount,  maxLength, FrequencyAdd, PhaseAdd)
%% Load pulse shape from .mat file
%
% The .mat file must contain the following variables:
%   ComplexAmplitude
%     Complex amplitude. Only used if Amplitude and Phase is not used (see
%     below).
%   Amplitude
%     Absolute value of pulse amplitude segments. Takes precedence over
%     ComplexAmplitude.
%   Phase
%     Phase of pulse segments in degrees. Takes precedence over
%     ComplexAmplitude.
%   nPulseSegments
%     Number of pulse segments. (Default: [])
%   DurationOfPulse
%     Total duration of shaped pulse in seconds. If not set, it is sum of
%     DurationOfPulseSegments (see below).
%   DurationOfPulseSegments
%     Duration of pulse segments in seconds. If not set, Segments have equal
%     duration.
%   CenterOfPulses
%     Effective center of NMR excitation in seconds. (Default: 0)
%   TimeOffsetToCenterOfPulse
%     Time from center of the pulse to the center of the effective NMR
%     excitation in seconds. Similar to CenterOfPulses. (Default: 0)
%   Bandwidth
%     Pulse bandwidth in Hz. (Default: 0)
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------


%% load file
load DemoPulseData.mat

%% check content of loaded file
if ~exist('ComplexAmplitude', 'var') || isempty(ComplexAmplitude) && ...
    exist('Amplitude', 'var') && ~isempty(Amplitude) && ...
    exist('Phase', 'var') && ~isempty(Phase)
  ComplexAmplitude = Amplitude .* exp(1i*Phase/180*pi);
end
if ~exist('Amplitude', 'var') || isempty(Amplitude)
  Amplitude = abs(ComplexAmplitude);
  Phase = angle(ComplexAmplitude) / pi * 180;
end
if ~exist('Phase', 'var') || isempty(Phase)
  Phase = zeros(size(Amplitude));
end
if ~exist('nPulseSegments', 'var') || isempty(nPulseSegments)
  nPulseSegments = numel(Amplitude);
end
if ~exist('FrequencyOffset', 'var') || isempty(FrequencyOffset)
  FrequencyOffset=zeros(size(Amplitude));
end
if (~exist('DurationOfPulse', 'var') || isempty(DurationOfPulse)) && ...
    exist('DurationOfPulseSegments', 'var') && ~isempty(DurationOfPulseSegments)
  DurationOfPulse = sum(DurationOfPulseSegments);
end
if ~exist('DurationOfPulseSegments', 'var') || isempty(DurationOfPulseSegments)
  DurationOfPulseSegments = DurationOfPulse/nPulseSegments;
end
if ~exist('DurationOfPulse', 'var') || isempty(DurationOfPulse)
  DurationOfPulse = sum(DurationOfPulseSegments);
end
if ~exist('CenterOfPulses', 'var') || isempty(CenterOfPulses)
  CenterOfPulses = 0;
end
if ~exist('TimeOffsetToCenterOfPulse', 'var') || isempty(TimeOffsetToCenterOfPulse)
  TimeOffsetToCenterOfPulse = 0;
end
if ~exist('Bandwidth', 'var') || isempty(Bandwidth)
  Bandwidth = 0;
end

if isscalar(DurationOfPulseSegments)
  DurationOfPulseSegments = ones(size(Amplitude)) * DurationOfPulseSegments;
end
if isscalar(Phase)
  Phase = ones(size(Amplitude)) * Phase;
end
if isscalar(FrequencyOffset)
  FrequencyOffset = ones(size(Amplitude)) * FrequencyOffset;
end

if nargin == 2
  % short path (additional syntax)
  if strcmp(Center, 'Amp')
    % factor in amp of the pulse to get a 90 degrees flip in the same time
    % as a 90 degrees rect pulse
    TXShape = NaN;              % PulseLength/HW.tFlip90Def
  elseif strcmp(Center, 'Time')
    % factor in time of the pulse to get the same BW as an rect pulse
    % (BW of rect pulse = 1/length)
    TXShape = NaN;              % PulseLength * BW
  elseif strcmp(Center, 'Bandwidth')
    TXShape = Bandwidth;        % BW, use to fix BW of pulse (set 'Time' to nan)
  elseif strcmp(Center, 'Length')
    TXShape = DurationOfPulse;  % Length, use to fix length of pulse  (set 'Amp' and 'Time' to nan)
  else
    TXShape = NaN;
  end
  return;
end

% Set bandwidth and length of your pulse for the slice calculation
% put Shape values into TXShape
TXShape(1).Start = CenterOfPulses - DurationOfPulse/2 - TimeOffsetToCenterOfPulse + ...
  cumsum(DurationOfPulseSegments(:)) - DurationOfPulseSegments(1) + Center;
TXShape(1).Center = CenterOfPulses + Center;
TXShape(1).CenterOffset = TimeOffsetToCenterOfPulse;
TXShape(1).Duration = DurationOfPulseSegments(:);
TXShape(1).Frequency = FrequencyOffset(:) + FrequencyAdd;
TXShape(1).Phase = Phase(:) + PhaseAdd;
% Use gamma that better matches the frequency of the pulse
% FIXME: Could this be an issue with (very) off-center slice pulses?
if abs(Pulse.Frequency - HW.fLarmorX) < abs(Pulse.Frequency - HW.fLarmor)
  TXShape(1).Amplitude = Amplitude(:) * 2*pi/HW.GammaX;
else
  TXShape(1).Amplitude = Amplitude(:) * 2*pi/HW.GammaDef;
end


if 0
  %% generate a mat file
  nPulseSegments = 100;                                         % Number of pulse segments
  DurationOfPulse = 200e-6;                                     % Duration of shaped pulse
  DurationOfPulseSegments = DurationOfPulse/nPulseSegments;     % Duration of pulse segment
  ComplexAmplitude = sinc(linspace(-2,2,nPulseSegments))*6000;  % complex Amplitude
  Amplitude = abs(ComplexAmplitude);                            % abs of pulse Amplitude
  Phase = angle(ComplexAmplitude)/pi*180;                       % phase of pulse in deg
  FrequencyOffset = 0;                                          % Frequency offset of pulse in deg
  CenterOfPulses = 0;                                           % Center of NMR excitation
  TimeOffsetToCenterOfPulse = 0;                                % Time from duration center to the center of NMR excitation
  Bandwidth = 10e3;                                             % pulse bandwidth. Set to zero if not used
  save('DemoPulseData.mat',...
       'nPulseSegments',...
       'DurationOfPulse',...
       'DurationOfPulseSegments',...
       'ComplexAmplitude',...
       'Amplitude',...
       'Phase',...
       'FrequencyOffset',...
       'CenterOfPulses',...
       'TimeOffsetToCenterOfPulse',...
       'Bandwidth');
end

if 0
  %% translate a PulseData.mat file
  load PulseData.mat                  % load pulse data with amp (Amplitude) and pha (Phase)
  Amplitude = amp;                    % Amplitude
  Phase = pha + 90;                   % Phase in degrees

  DurationOfPulse = 0.5e-3;           % Set duration of the shaped pulse

  CenterOfPulses = 0;                 % Center of NMR excitation
  TimeOffsetToCenterOfPulse = 240e-6; % Time from duration center to the center of NMR excitation

  Bandwidth = 20e3;                   % pulse bandwidth. Set to zero if not used

  save('PulseData2.mat', ...
       'DurationOfPulse', ...
       'Amplitude', ...
       'Bandwidth', ...
       'CenterOfPulses', ...
       'TimeOffsetToCenterOfPulse', ...
       'Phase');
end

end
