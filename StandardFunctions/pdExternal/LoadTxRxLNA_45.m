%% Settings for LNA broadband experimental 

if ~exist('iDevice', 'var'), iDevice = 1; end

HW.RX.LnaSN = 45;

HW.TX(iDevice).Uout2PaUout(2) = HW.TX(iDevice).Uout2PaUout(2) * 10^(-0.2/20);  % -0.2 dB Gain (TX to Coil) @ 24.3 MHz
HW.TX(iDevice).ChannelDef = 2;  % default TX rf channel
HW.RX2TXdeadTime = 10e-6;  % dead-time of receiver before TX pulse
HW.TX2RXdeadTime = 50e-6;  % dead-time of receiver after TX pulse
HW.TX(iDevice).BlankOffset = 9e-6;  % >7 us - blank of transmit before TX pulse
HW.TX(iDevice).BlankPostset = 5e-6;  % blank of transmit after TX pulse

HW.TX(iDevice).BlankOffsetAQ = 2e-6;  % blank of receiver before TX pulse
HW.TX(iDevice).BlankPostsetAQ = 5e-6;  % blank of receiver after TX pulse
HW.TX(iDevice).BlankAQ = 1;  % switch TRx to 50 Ohm resistor during TX pulse, to avoid saturation.

HW.RX.VGAGainDef = HW.RX.VGAGainMax/3;  % default receiver gain

HW.RX.LNAGain = 10^(((-51.3934)-(-71.658))/20);  % 20.2646 dB gain @ 13100261.232 MHz F=0.80054 dB (-71.658 dBm cal) with 13.1 MHz passiv switch 

HW.TX(iDevice).Max.PaUout(2) = min(HW.TX(iDevice).Max.PaUout(2), 10);  % max transmit voltage
HW.TX(iDevice).Def.PaUout(2) = min(HW.TX(iDevice).Def.PaUout(2), 4);  % def transmit voltage

LoadExtLNA_Cal;
