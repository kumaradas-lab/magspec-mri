function isMatching = checkDeviceSerial(HW, deviceSerial, caller)
%% Check if the connected console matches a serial number
%
%   isMatching = checkDeviceSerial(HW, deviceSerial, caller)
%
% INPUT:
%
%   HW
%       HW object or structure.
%
%   deviceSerial
%       Device serial number of the console that is expected to be currently
%       connected.
%
%   caller
%       String which identifies from where this function is called.
%       (Default: '')
%
%
% OUTPUT:
%
%   isMatching
%       Boolean value that indicates whether the currently connected device has
%       the expected serial number (true) or not (false).
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% check input
if nargin < 2
  error('PD:checkDeviceSerial:InvalidArgs', ...
    'Must be called with at least two input arguments.');
end

if ~isa(HW, 'PD.HWClass') ...
    && (~isstruct(HW) || isemptyfield(HW, {'MMRT', 'DeviceSerial'}))
  error('PD:checkDeviceSerial:InvalidArgs', ...
    'Invalid input argument "HW".');
end

if ~isnumeric(deviceSerial)
  error('PD:checkDeviceSerial:InvalidArgs', ...
    'Invalid input argument "deviceSerial".');
end


%% default parameters
if nargin < 3
  caller = '';
end

if ~ischar(caller)
  error('PD:checkDeviceSerial:InvalidArgs', ...
    'Invalid input argument "caller".');
end


%% check device serial
% FIXME: Handle multiple connected devices?
isMatching = (HW.MMRT(1).DeviceSerial == deviceSerial);

if ~isMatching
  if isempty(caller)
    callerStr = '';
  else
    callerStr = sprintf(' or re-create the file "%s"', caller);
  end
  warningCommand = 'warning(''off'', ''PD:checkDeviceSerial:NoMatch'');';
  warning('PD:checkDeviceSerial:NoMatch', ...
    ['The expected device serial "%d" doesn''t match the serial number of the connected device "%d".\n', ...
    'Check if the correct device is connected%s.\n', ...
    'Alternatively, set <a href="matlab: %s disp(''%s'')">%s</a> to suppress this warning.\n'], ...
    deviceSerial, HW.MMRT(1).DeviceSerial, callerStr, ...
    warningCommand, regexprep(warningCommand, '''', ''''''), warningCommand);
end

end
