function gain = get_CalibratedGain(HW, Network)
%% Correct gain using cable calibration and latency
%
%   gain = get_CalibratedGain(HW, Network)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------


CableGainRaw = interp1(HW.NetworkCal.Cable.FrequencyGain, ...
  HW.NetworkCal.Cable.GainRaw, Network.FrequencyGain);

AddCableLatency = HW.AddCableLength / (HW.Constant.SpeedOfLight*HW.RelativCableSpeed) * 2;
gain = abs(Network.GainRaw) ./ abs(CableGainRaw) .* ...
  exp(-1i*((angle(Network.GainRaw)-angle(CableGainRaw)) - 2*pi*AddCableLatency*Network.FrequencyGain));

% gain = abs(Network.GainRaw) .* ...
%   exp(-1i*((angle(Network.GainRaw)) - 2*pi*AddCableLatency*Network.FrequencyGain));

% Reflection = (real(r)-1i.*imag(r)) .* exp(+1i*2*pi*AddCableLatency*Network.Frequency);
% Reflection = r .* exp(+1i*2*pi*AddCableLatency*Network.Frequency);
% Reflection = (r.')' .* exp(+1i*2*pi*AddCableLatency*Network.Frequency);

end
