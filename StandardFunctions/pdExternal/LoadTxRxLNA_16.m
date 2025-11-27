%% Settings for LNA S/N 16

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 16;

% LNA only (no switch)
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel

HW.RX(iDevice).VGAGainDef = HW.RX.VGAGainMax/3;  % reduce VGA gain to avoid saturation
HW.RX(iDevice).LNAGain = 10^(((-58.45)-(-80))/20);  % 21.55 dB gain @ 24.0 MHz F=0.9 dB

% %% TRx switch during transmit at Tx2
% HW.TX(iDevice).BlankOffsetAQ = 1000e-9;  % blank of receiver before TX pulse
% HW.TX(iDevice).BlankPostsetAQ = 4000e-9;  % blank of receiver after TX pulse
% HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.
