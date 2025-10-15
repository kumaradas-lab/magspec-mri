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
%           At least one of the following input properties must be set:
%           'Norm', 'Mmrt', 'Uout', 'PaUout', 'Amplitude'
%           If more than one of the former input properties is set, the
%           strictest applies.
%           Additionally, the following properties can be set. If they are
%           omitted, default values are used:
%           'Frequency' (default: HW.fLarmor)
%           'Device'    (default: 1)
%           'Channel'   (default: HW.TX(Device).ChannelDef)
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
% (C) Copyright 2020-2024 Pure Devices GmbH, Wuerzburg, Germany
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

% parse varargin
inputProps = {'Norm', 'Mmrt', 'Uout', 'PaUout', 'Amplitude'};
frequency = [];
iDevice = [];
txChannel = [];
inputPropIdx = 0;
for iArgin = 1:2:numel(varargin)
  if ~ischar(varargin{iArgin})
    error('PD:get_TX_Amplitude:InvalidProp', ...
      'Properties must be a string.');
  end
  switch lower(varargin{iArgin})
    case lower(inputProps)
      inputPropIdx = inputPropIdx + 1;
      inputProp{inputPropIdx} = varargin{iArgin};  %#ok<AGROW>
      inputVal{inputPropIdx} = varargin{iArgin+1};  %#ok<AGROW>
    case 'frequency'
      frequency = varargin{iArgin+1};
    case 'device'
      iDevice = varargin{iArgin+1};
    case 'channel'
      txChannel = varargin{iArgin+1};
    otherwise
      error('PD:get_TX_Amplitude:UnknownProp', ...
        'Unknown input property "%s".', varargin{iArgin});
  end
end

% default input
if isempty(frequency)
  frequency = HW.fLarmor;
end
if isempty(iDevice)
  iDevice = 1;
end
if isempty(txChannel)
  txChannel = HW.TX(iDevice).ChannelDef;
end

if inputPropIdx == 0
  error('PD:get_TX_Amplitude:NoInputVal', ...
    'One input value must be set as property value pair.');
end


%% use helper object to get corresponding amplitude
% get correct lower/upper case of properties
lowerInputProps = cellfun(@lower, inputProps, 'UniformOutput', false);
lowerInputProp = cellfun(@lower, inputProp, 'UniformOutput', false);
[~, iLetterCase, iInput] = intersect(lowerInputProps, lowerInputProp);
idxSorted = 1:numel(inputProp);
idxSorted(iInput) = idxSorted;
inputProp = inputProps(iLetterCase(idxSorted));

ampOk = true;
% get maximum size of all input values
inputValSize = cellfun(@size, inputVal, 'UniformOutput', false);
inputValSize = max(cat(3, inputValSize{:}), [], 3);
inputVal = cellfun(@(x) x .* ones(inputValSize), inputVal, 'UniformOutput', false);
if isscalar(frequency)
  % scalar frequency allows creating the helper object only once
  TXCal = PD.TXMaxDef(HW.TX(iDevice), 'Cal', frequency);
  amplitude = zeros(inputValSize);

  for iVal = 1:numel(amplitude)
    for iProp = 1:inputPropIdx
      % FIXME: Is a loop good enough?
      if isnan(inputVal{iProp}(iVal))
        amplitude(iVal) = NaN;
        continue;
      end

      TXCal.(inputProp{iProp})(txChannel) = inputVal{iProp}(iVal);
    end

    if (inputPropIdx==1) && ampOk ...
        && TXCal.([inputProp{iProp}, 'Calibrated'])(txChannel) < 0.999 * inputVal{iProp}(iVal)
      % FIXME: Should we check if all(!) input properties "failed" in case of
      % more than one input property?
      ampOk = false;
    end

    amplitude(iVal) = TXCal.AmplitudeCalibrated(txChannel);
  end

else
  % non-scalar frequency: create helper object for each loop iteration

  % check common size
  commonSize = ones(max(inputValSize, size(frequency)));
  inputVal = cellfun(@(x) x .* commonSize, inputVal, 'UniformOutput', false);
  frequency = frequency .* commonSize;
  amplitude = zeros(size(commonSize));

  % FIXME: Is a loop good enough?
  for iVal = 1:numel(amplitude)
    TXCal = PD.TXMaxDef(HW.TX(iDevice), 'Cal', frequency(iVal));

    for iProp = 1:inputPropIdx
      if isnan(inputVal{iProp}(iVal))
        amplitude(iVal) = NaN;
        continue;
      end

      TXCal.(inputProp{iProp})(txChannel) = inputVal{iProp}(iVal);
    end

    if (inputPropIdx==1) && ampOk ...
        && (TXCal.([inputProp{iProp}, 'Calibrated'])(txChannel) < 0.999 * inputVal{iProp}(iVal))
      % FIXME: Should we check if all(!) input properties "failed" in case of
      % more than one input property?
      ampOk = false;
    end

    amplitude(iVal) = TXCal.AmplitudeCalibrated(txChannel);
  end
end


%% check result
isZeroInputVal = cellfun(@eq, inputVal, repmat({0}, 1, inputPropIdx), 'UniformOutput', false);
if any(~any(cat(3, isZeroInputVal{:}), 3) & (amplitude == 0)) || ~(ampOk || inputPropIdx>1)
  warning('PD:get_TX_Amplitude:InvalidInput', ...
    ['Could not calculate amplitude for (some of) the selected input. ', ...
    'Are the requested parameters outside the calibrated range?']);
end


end
