%% Settings for LNA S/N 11

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 11;

% Resonant coil 15 mm
HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.2/20);
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 5e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 30e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 160e-9;  % blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 160e-9;  % blank of transmit after TX pulse

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % default receiver gain

HW.RX(iDevice).LNAGain = 10^(22.38/20); % 22.38 dB gain @ 24.71 MHz F=1 dB

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 100);

%% TRx switch during transmit at Tx2
HW.TX(iDevice).BlankOffsetAQ = 800e-9;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 2000e-9;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % wwitch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.
