function deadtime = get_DeadTimeRX2TX(HW, fSample, iDevice)
%% Calculate dead time to the next TX pulse after an acquisition
%
%   deadtime = get_DeadTimeRX2TX(HW, fSample)
%
% Determine reasonable dead time before TX pulse taking the properties of the
% down-sampler at the given sampling rate into acount.
%
%
% INPUT:
%
%   HW
%       HW object or structure.
%
%   fSample
%       Sampling frequency in Hz of the acquisition window preceeding the rf
%       pulse.
%
%
% OUTPUT:
%
%   deadtime
%       Minimum time in seconds required between the end of an acquisition
%       window with the given sampling frequency and the start of a following rf
%       pulse.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------


%% default input
if nargin < 3, iDevice = 1; end


%% minimum dead time of receiver before rf pulse
if numel(HW.TX2RXdeadTime) >= iDevice
  RX2TXdeadTime = HW.RX2TXdeadTime(iDevice);
else
  RX2TXdeadTime = HW.RX2TXdeadTime(1);
end


%% potentially increase dead time for blanking the AQ before the rf pulse
% FIXME: Do we need to take the SamplingFactor into account here?
% This is probably good enough as long as the SamplingFactor is 1 or larger than
% approximately 4.
deadtime = max(RX2TXdeadTime, HW.TX(iDevice).BlankOffset) + 4 ./ fSample;

end
