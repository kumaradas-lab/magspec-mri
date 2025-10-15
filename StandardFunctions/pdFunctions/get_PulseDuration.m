function pulseDuration = get_PulseDuration(HW, flipAngle, amplitude)
%% Calculate RF pulse duration using RF amplitude and target flip angle
%
%     pulseDuration = get_PulseDuration(HW, flipAngle, amplitude)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

iDevice = 1;  % FIXME: Support multiple MMRT devices
if nargin == 2
  amplitude = HW.TX(iDevice).AmpDef;
end

pulseDuration = 1 ./ (amplitude .* HW.TX(iDevice).Amp2Hz) ./ (360) .* flipAngle;

end
