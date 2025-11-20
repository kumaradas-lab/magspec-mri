%% Settings for LNA #24 24.377452 MHz

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 24;

% HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.4/20);  % -0.4 dB Gain @ 23.5 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 1e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 4e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 1e-6;  % blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 480e-9;  % blank of transmit after TX pulse
HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/6;  % default receiver gain

HW.TX(iDevice).BlankOffsetAQ = 400e-9;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 400e-9;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX(iDevice).LNAGain = 10^(((-36.7699)-(-60))/20);  % 23.2301 dB gain @ 24377838.756 MHz F=1.0667 dB 50 Ohm (-60 dBm cal)

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 100);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 50);  % def transmit voltage
