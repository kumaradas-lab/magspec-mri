function deadtime = get_DeadTimeTX2RX(HW, fSample, iDevice)
%% Calculate dead time of receiver after TX pulse
%
%   deadtime = get_DeadTimeTX2RX(HW, fSample)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

if nargin < 3, iDevice = 1; end

if HW.TX(iDevice).DampCoil.Enable
  TX2RXdeadTime = HW.TX(iDevice).DampCoil.TX2RXdeadTime;
else
  if numel(HW.TX2RXdeadTime) >= iDevice
    TX2RXdeadTime = HW.TX2RXdeadTime(iDevice);
  else
    TX2RXdeadTime = HW.TX2RXdeadTime(1);
  end
end

% FIXME: Do we need to take the SamplingFactor into account here?
% This is probably good enough as long as the SamplingFactor is 1 or larger than
% approximately 4.
deadtime = max(TX2RXdeadTime, HW.TX(iDevice).BlankPostset) + 4/fSample;

end
