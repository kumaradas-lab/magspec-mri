function amp = get_RX_Amplitude(HW, amp, outType, varargin)
%% Query value of RX amplitude at different points in amplification chain
%
%    ampOut = get_RX_Amplitude(HW, ampIn, outType, 'prop', val, ...)
%
%
% INPUT:
%
%   HW
%           HW object
%
%   ampIn
%           Received signal amplitude (e.g. as in data.data). By default, these
%           amplitudes are measured in T (see settings in HW.RX).
%
%   outType
%           Identifier for position in the amplification chain. The following
%           types are supported:
%           'LnaUin', 'Uin', 'VgaUin', 'VgaUout', 'MmrtUin', 'AdcUin',
%           'Norm', 'Adc'
%
%   'prop', val, ...
%           List of property value pairs.
%           The following properties can be set additionally. If they are
%           omitted, default values are used:
%           'Device'    (default: 1)
%           'AQ'        (default: AQ.Gain = HW.RX.GainDef;
%                                 AQ.Frequency = HW.fLarmor;)
%
%
% OUTPUT:
%
%   ampOut
%           The amplitude at the selected position in the amplification chain.
%
%
% EXAMPLE:
%
%   Demo_FID
%   dataUin = get_RX_Amplitude(SeqOut.HW, data.data, 'Uin', 'AQ', SeqOut.AQ);
%
% ------------------------------------------------------------------------------
% (C) Copyright 2020-2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% check and parse input
if nargin < 3
  error('PD:query_RX_Amplitude:MissingArgs', ...
    'The first three input arguments must be set.');
end

if ~isa(HW, 'PD.HWClass') && ~isstruct(HW)
  error('PD:query_RX_Amplitude:InvalidHW', ...
    'The first input argument must be a PD.HWClass object or struct.');
end

if ~isnumeric(amp)
  error('PD:query_RX_Amplitude:InvalidAmpIn', ...
    '"ampIn" must be a numeric array.');
end

if ~ischar(outType)
  error('PD:query_RX_Amplitude:InvalidOutType', ...
    '"outType" must be a string.');
end

ampChain = {'LnaUin', 'Uin', 'VgaUin', 'VgaUout', 'MmrtUin', 'AdcUin', 'Norm', 'Adc'};
outPos = find(strcmp(outType, ampChain));

if isempty(outPos)
  error('PD:query_RX_Amplitude:InvalidOutType', ...
    '"outType" must identify a position in the amplification chain.');
end

% get optional input
if mod(numel(varargin), 2)
  error('PD:query_RX_Amplitude:InvalidPropVal', ...
    'Additional arguments must be property value pairs.');
end

%% default values
iDevice = 1;
AQ = struct();

for iArgin = 1:2:numel(varargin)
  if ~ischar(varargin{iArgin})
    error('PD:query_RX_Amplitude:InvalidProp', ...
      'Properties must be a string.');
  end
  switch lower(varargin{iArgin})
    case 'aq'
      AQ = varargin{iArgin+1};
    case 'device'
      iDevice = varargin{iArgin+1};
    otherwise
      error('PD:query_RX_Amplitude:UnknownProp', ...
        'Unknown input property "%s".', varargin{iArgin});
  end
end


%% convert to selected position in the amplification chain and return
% LnaUin
amp = amp * HW.RX(iDevice).Amplitude2LnaUin;
if outPos < 2
  return;
end

% Uin
amp = amp * HW.RX(iDevice).LNAGain;
if outPos < 3
  return;
end

% VgaUin
amp = amp * HW.RX(iDevice).Uin2VgaUin;
if outPos < 4
  return;
end

% VgaUout
gain = HW.RX(iDevice).GainDef;
if isfield(AQ, 'Gain')
  gain = AQ.Gain;
end

vgaGain = gain ./ ...
  (HW.RX(iDevice).Amplitude2LnaUin .* HW.RX(iDevice).LNAGain .* ...
  HW.RX(iDevice).Uin2VgaUin .* HW.RX(iDevice).VgaUout2MmrtUin .* ...
  HW.RX(iDevice).MmrtUin2AdcUin .* HW.RX(iDevice).AdcUin2Norm);

amplitudeCalGain = 1;
if (isa(HW.RX(iDevice), 'PD.RX') || isfield(HW.RX(iDevice), 'Calibration')) && ...
    ~isempty(HW.RX(iDevice).Calibration) && ...
    isfield(HW.RX(iDevice).Calibration, 'Frequency') && ...
    ~isempty(HW.RX(iDevice).Calibration.Frequency)
  frequency = HW.fLarmor;
  if isfield(AQ, 'Frequency')
    frequency = AQ.Frequency;
  end

  amplitudeCalGain = interp1(HW.RX(iDevice).Calibration.Frequency, ...
    HW.RX(iDevice).Calibration.Gain, frequency, 'pchip', 1);
end

amp = bsxfun(@times, amp * vgaGain, amplitudeCalGain);
if outPos < 5
  return;
end

% MmrtUin
amp = amp * HW.RX(iDevice).VgaUout2MmrtUin;
if outPos < 6
  return;
end

% AdcUin
amp = amp * HW.RX(iDevice).MmrtUin2AdcUin;
if outPos < 7
  return;
end

% Norm
amp = amp * HW.RX(iDevice).AdcUin2Norm;
if outPos < 8
  return;
end

% Adc
amp = amp * HW.RX(iDevice).Norm2Adc;

end
