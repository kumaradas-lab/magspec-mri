%% Settings for LNA #48 24.3 MHz

% Low power broadband LNA built Aug 2023

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX(iDevice).LnaSN = 48;

HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.2/20);  % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 10e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 9e-6;  % >7 us - blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 5e-6;  % blank of transmit after TX pulse

HW.TX(iDevice).BlankOffsetAQ = 2e-6;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 5e-6;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX(iDevice).VGAGainDef = HW.RX(iDevice).VGAGainMax/3;  % default receiver gain

HW.RX(iDevice).LNAGain = 10^(((-49.4223)-(-71.1031))/20);  % 21.6808 dB gain @ 24300454.0347 MHz F=1.1289 dB (-71.1031 dBm cal)
% HW.RX(iDevice).LNAGain = 10^(((-49.0633)-(-70.7971))/20);  % 21.7338 dB gain @ 21490393.8015 MHz F=1.3551 dB (-70.7971 dBm cal)

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 10);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 4);  % def transmit voltage

LoadExtLNA_Cal;
