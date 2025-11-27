function deadtime = get_DeadTimeTX2RX(HW, fSample, iDevice)
%% Calculate dead time of receiver after TX pulse
%
%   deadtime = get_DeadTimeTX2RX(HW, fSample)
%
% Determine reasonable dead time after TX pulse taking the properties of the
% down-sampler at the given sampling rate into acount.
%
%
% INPUT:
%
%   HW
%       HW object or structure.
%
%   fSample
%       Sampling frequency in Hz of the acquisition window following the rf
%       pulse.
%
%
% OUTPUT:
%
%   deadtime
%       Minimum time in seconds required between the end of an rf pulse and the
%       start of a following acquisition window with the given sampling
%       frequency.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------


%% default input
if nargin < 3, iDevice = 1; end


%% minimum dead time for ring down (with or without active coil damping)
if HW.TX(iDevice).DampCoil.Enable
  TX2RXdeadTime = HW.TX(iDevice).DampCoil.TX2RXdeadTime;
else
  if numel(HW.TX2RXdeadTime) >= iDevice
    TX2RXdeadTime = HW.TX2RXdeadTime(iDevice);
  else
    TX2RXdeadTime = HW.TX2RXdeadTime(1);
  end
end


%% potentially increase dead time for unblanking the AQ after the rf pulse
% FIXME: Do we need to take the SamplingFactor into account here?
% This is probably good enough as long as the SamplingFactor is 1 or larger than
% approximately 4.
deadtime = max(TX2RXdeadTime, HW.TX(iDevice).BlankPostset) + 4 ./ fSample;

end
