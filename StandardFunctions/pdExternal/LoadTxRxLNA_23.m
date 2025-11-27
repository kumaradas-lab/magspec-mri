%% Settings for LNA #23

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 23;

HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.4/20);  % -0.4 dB Gain @ 23.5 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 3e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 13e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 2e-6;  % blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 480e-9;  % blank of transmit after TX pulse
HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % default receiver gain

HW.TX(iDevice).BlankOffsetAQ = 2000e-9;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 4000e-9;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX(iDevice).LNAGain = 10^(21.3/20);  % 21.3 dB LNA receive gain @ 23.5 MHz

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 100);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 100);  % def transmit voltage
