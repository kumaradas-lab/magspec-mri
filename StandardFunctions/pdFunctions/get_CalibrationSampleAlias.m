function alias = get_CalibrationSampleAlias(HW, uid, iDevice)
%% Get the alias (display name) for calibration sample from its UID
%
%   alias = get_CalibrationSampleAlias(HW, uid, iDevice)
%
%
% INPUT:
%
%   HW
%       HW object or structure.
%
%   uid
%       String or cell array of strings with UIDs of calibration samples for
%       which data is available in the file
%       HW.RX(iDevice).LoadCalibrationSamplePath.
%
%   iDevice
%       Index in case multiple devices are connected.
%       (Default: 1)
%
%
% OUTPUT:
%
%   alias
%       Corresponding alias for the input uid. If the input uid is a string,
%       alias is also a string. If the input uid is a cell string, alias is also
%       a cell string.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input parameters
if nargin < 3 || isempty(iDevice)
  iDevice = 1;
end


%% load calibration sample data from file
CalibrationSample = LoadCalibrationSample(HW, iDevice);


%% look up alias for each UID
if iscell(uid)
  alias = cell(size(uid));
  for iCell = 1:numel(uid)
    sample = CalibrationSample(strcmp(CalibrationSample.UID, uid{iCell}), :).Alias;
    if isempty(sample)
      alias(iCell) = {[uid{iCell}, ' (unknown)']};
    else
      alias(iCell) = sample(1,1);
    end
  end
else
  sample = CalibrationSample(strcmp(CalibrationSample.UID, uid), :).Alias;
  if isempty(sample)
    alias = 'unknown';
  else
    alias = sample(1,1);
  end
end

end
