function deadtime = get_DeadTimeRX2TX(HW, fSample, iDevice)
%% Calculate dead time to the next pulse the receiver needs after an acquisition
%
%   deadtime = get_DeadTimeRX2TX(HW, fSample)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

if nargin < 3, iDevice = 1; end

if numel(HW.TX2RXdeadTime) >= iDevice
  RX2TXdeadTime = HW.RX2TXdeadTime(iDevice);
else
  RX2TXdeadTime = HW.RX2TXdeadTime(1);
end

% FIXME: Do we need to take the SamplingFactor into account here?
% This is probably good enough as long as the SamplingFactor is 1 or larger than
% approximately 4.
deadtime = max(RX2TXdeadTime, HW.TX(iDevice).BlankOffset) + 4 ./ fSample;

end
