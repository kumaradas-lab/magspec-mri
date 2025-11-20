%% Settings for LNA #25 22.813840  MHz

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 25;

% HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2)*10^(-0.4/20);  % -0.4 dB Gain @ 23.5 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 1e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 5e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 1e-6;  % blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 1.4e-6;  % blank of transmit after TX pulse
HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/6;  % default receiver gain

HW.TX(iDevice).BlankOffsetAQ = 1.4e-6;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 1e-6;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

% HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/6;  % default receiver gain

HW.RX(iDevice).LNAGain = 10^(((-35.599)-(-60))/20);  % 24.401 dB gain @ 22814173.6183 MHz F=0.85267 dB (-60 dBm cal)

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 101);  % max transmit voltage
% HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 50);  % def transmit voltage
