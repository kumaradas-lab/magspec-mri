%% Settings for LNA S/N 15

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 15;

% Resonant coil 15 mm
HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.2/20);
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 5e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 40e-6;  % 40 us (slow PIN diodes in LNA #12) - dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 160e-9;  % blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 160e-9;  % blank of transmit after TX pulse

if 0
  HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % default receiver gain
  HW.RX(iDevice).LNAGain = 10^(((-54.4)-(-80))/20);  % 26.0 dB gain @ 13.0 MHz F=1 dB
else
  HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/1;  % default receiver gain
  HW.RX(iDevice).LNAGain = 10^(((-54.28)-(-80))/20);  % 26.0 dB gain @ 13.0 MHz F=0.6 dB
end

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 100);

%% TRx switch during transmit at Tx2
HW.TX(iDevice).BlankOffsetAQ = 800e-9;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 2000e-9;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.
