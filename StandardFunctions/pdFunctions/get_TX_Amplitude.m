function amplitude = get_TX_Amplitude(HW, varargin)
%% Get TX amplitude that corresponds to a given input value
%
%   amplitude = get_TX_Amplitude(HW, 'prop', val, ...)
%
%
% INPUT:
%
%   HW
%           HW object
%
%   'prop', val, ...
%           List of property value pairs.
%           Exactly one of the following input properties must be set:
%           'Norm', 'Mmrt', 'Uout', 'PaUout', 'Amplitude'
%           Additionally, the following properties can be set. If they are
%           omitted, default values are used:
%           'Frequency' (default: HW.fLarmor)
%           'Channel'   (default: HW.TX.ChannelDef)
%
%
% OUTPUT:
%
%   amplitude
%           The rf amplitude in Tesla that corresponds to the input values.
%           Calibrations (factors, maps, ...) are taken into account.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% check and parse input
if mod(numel(varargin), 2)
  error('PD:get_TX_Amplitude:InvalidPropVal', ...
    'Additional arguments must be property value pairs.');
end

if numel(varargin) < 2
  error('PD:get_TX_Amplitude:NoInputVal', ...
    'One input value must be set as property value pair.');
end

% default input
frequency = HW.fLarmor;
txChannel = HW.TX.ChannelDef;

% parse varargin
inputProps = {'Norm', 'Mmrt', 'Uout', 'PaUout', 'Amplitude'};
foundInputVal = false;
for iArgin = 1:2:numel(varargin)
  if ~ischar(varargin{iArgin})
    error('PD:get_TX_Amplitude:InvalidProp', ...
      'Properties must be a string.');
  end
  switch lower(varargin{iArgin})
    case lower(inputProps)
      if foundInputVal
        error('PD:get_TX_Amplitude:MultipleInputVal', ...
          'Only one input value can be set at a time.');
      end
      inputProp = varargin{iArgin};
      inputVal = varargin{iArgin+1};
      foundInputVal = true;
    case 'frequency'
      frequency = varargin{iArgin+1};
    case 'channel'
      txChannel = varargin{iArgin+1};
    otherwise
      error('PD:get_TX_Amplitude:UnknownProp', ...
        'Unknown input property "%s".', varargin{iArgin});
  end
end

if ~foundInputVal
  error('PD:get_TX_Amplitude:NoInputVal', ...
    'One input value must be set as property value pair.');
end

%% create helper object
TXCal = PD.TXMaxDef(HW.TX, 'Cal', frequency);

% get correct lower/upper case of properties
inputProp = inputProps{strcmpi(inputProps, inputProp)};

% FIXME: Is a loop good enough?
amplitude = zeros(size(inputVal));
for iVal = 1:numel(inputVal)
  if isnan(inputVal(iVal))
    amplitude(iVal) = NaN;
    continue;
  end

  TXCal.(inputProp)(txChannel) = inputVal(iVal);

  amplitude(iVal) = TXCal.AmplitudeCalibrated(txChannel);
end

end
